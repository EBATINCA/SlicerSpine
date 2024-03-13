import os
import unittest
import logging
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin

#
# SegmentAxisAlignment
#

class SegmentAxisAlignment(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Segment Axis Alignment"
    self.parent.categories = ["Spine Strength Prediction"]
    self.parent.dependencies = ['Data', 'Segmentations']
    self.parent.contributors = ["Csaba Pinter (EBATINCA), Csaba Hegedus (Win95 Kft)"]
    self.parent.helpText = """
This module calculates the angles corresponding to the minimum bounding box of a segment thus aligning it with the principal axes.
"""
    self.parent.acknowledgementText = """
This file was originally developed by Csaba Pinter (EBATINCA) and was funded by Beth Israel Deaconess Medical Center, Boston, MA, USA.
"""

#
# SegmentAxisAlignmentWidget
#

class SegmentAxisAlignmentWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent=None):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.__init__(self, parent)
    VTKObservationMixin.__init__(self)  # needed for parameter node observation
    self.logic = None
    self._parameterNode = None
    self._updatingGUIFromParameterNode = False

  def setup(self):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.setup(self)

    # Load widget from .ui file (created by Qt Designer).
    # Additional widgets can be instantiated manually and added to self.layout.
    uiWidget = slicer.util.loadUI(self.resourcePath('UI/SegmentAxisAlignment.ui'))
    self.layout.addWidget(uiWidget)
    self.ui = slicer.util.childWidgetVariables(uiWidget)

    # Set scene in MRML widgets. Make sure that in Qt designer the top-level qMRMLWidget's
    # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
    # "setMRMLScene(vtkMRMLScene*)" slot.
    uiWidget.setMRMLScene(slicer.mrmlScene)

    # Create logic class. Logic implements all computations that should be possible to run
    # in batch mode, without a graphical user interface.
    self.logic = SegmentAxisAlignmentLogic()

    # Hide widgets that are not used at the moment
    self.ui.supersamplingFactorLabel.visible = False
    self.ui.supersamplingFactorSpinBox.visible = False
    self.ui.subsampleAfterRotationCheckBox.visible = False

    # Connections

    # These connections ensure that we update parameter node when scene is closed
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

    # These connections ensure that whenever user changes some settings on the GUI, that is saved in the MRML scene
    # (in the selected parameter node).
    self.ui.inputSegmentSelectorWidget.currentNodeChanged.connect(self.updateParameterNodeFromGUI)
    self.ui.inputSegmentSelectorWidget.currentSegmentChanged.connect(self.updateParameterNodeFromGUI)
    self.ui.inputSegmentSelectorWidget.currentNodeChanged.connect(self.onInputChanged)
    self.ui.inputSegmentSelectorWidget.currentSegmentChanged.connect(self.onInputChanged)
    self.ui.outputTransformSelectorComboBox.currentNodeChanged.connect(self.updateParameterNodeFromGUI)
    self.ui.outputSegmentationSelector.currentNodeChanged.connect(self.updateParameterNodeFromGUI)
    self.ui.checkBox_RotateAP.toggled.connect(self.updateParameterNodeFromGUI)
    self.ui.inputVolumeSelector.currentNodeChanged.connect(self.updateParameterNodeFromGUI)
    self.ui.subsampleAfterRotationCheckBox.toggled.connect(self.updateParameterNodeFromGUI)
    self.ui.outputVolumeSelector.currentNodeChanged.connect(self.updateParameterNodeFromGUI)

    # Buttons
    self.ui.calculateButton.clicked.connect(self.onCalculateButtonClicked)
    self.ui.alignVolumeButton.clicked.connect(self.onAlignVolumeButtonClicked)
    self.ui.feedInputToStrengthModuleButton.clicked.connect(self.onFeedInputToStrengthModuleButtonClicked)

    # Hide feed input to strength module button if the module is not present
    if not hasattr(slicer.modules, 'vertebralstrengthcalculation'):
      self.ui.feedInputToStrengthModuleButton.visible = False

    # Make sure parameter node is initialized (needed for module reload)
    self.initializeParameterNode()

  def cleanup(self):
    """
    Called when the application closes and the module widget is destroyed.
    """
    self.removeObservers()

  def enter(self):
    """
    Called each time the user opens this module.
    """
    # Make sure parameter node exists and observed
    self.initializeParameterNode()

  def exit(self):
    """
    Called each time the user opens a different module.
    """
    # Do not react to parameter node changes (GUI wlil be updated when the user enters into the module)
    self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

  def onSceneStartClose(self, caller, event):
    """
    Called just before the scene is closed.
    """
    # Parameter node will be reset, do not use it anymore
    self.setParameterNode(None)

  def onSceneEndClose(self, caller, event):
    """
    Called just after the scene is closed.
    """
    # If this module is shown while the scene is closed then recreate a new parameter node immediately
    if self.parent.isEntered:
      self.initializeParameterNode()

  def initializeParameterNode(self):
    """
    Ensure parameter node exists and observed.
    """
    # Parameter node stores all user choices in parameter values, node selections, etc.
    # so that when the scene is saved and reloaded, these settings are restored.

    self.setParameterNode(self.logic.getParameterNode())

    # Select default input nodes if nothing is selected yet to save a few clicks for the user
    # if not self._parameterNode.GetNodeReference("InputVolume"):
    #   firstVolumeNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLScalarVolumeNode")
    #   if firstVolumeNode:
    #     self._parameterNode.SetNodeReferenceID("InputVolume", firstVolumeNode.GetID())

  def setParameterNode(self, inputParameterNode):
    """
    Set and observe parameter node.
    Observation is needed because when the parameter node is changed then the GUI must be updated immediately.
    """
    if inputParameterNode:
      self.logic.setDefaultParameters(inputParameterNode)

    # Unobserve previously selected parameter node and add an observer to the newly selected.
    # Changes of parameter node are observed so that whenever parameters are changed by a script or any other module
    # those are reflected immediately in the GUI.
    if self._parameterNode is not None:
      self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)
    self._parameterNode = inputParameterNode
    if self._parameterNode is not None:
      self.addObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

    # Initial GUI update
    self.updateGUIFromParameterNode()

  def updateGUIFromParameterNode(self, caller=None, event=None):
    """
    This method is called whenever parameter node is changed.
    The module GUI is updated to show the current state of the parameter node.
    """

    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return

    # Make sure GUI changes do not call updateParameterNodeFromGUI (it could cause infinite loop)
    self._updatingGUIFromParameterNode = True

    # Update node selectors and sliders
    self.ui.inputSegmentSelectorWidget.setCurrentNode(self._parameterNode.GetNodeReference("InputSegmentationNode"))
    self.ui.inputSegmentSelectorWidget.setCurrentSegmentID(self._parameterNode.GetParameter("InputSegmentID"))
    self.ui.outputTransformSelectorComboBox.setCurrentNode(self._parameterNode.GetNodeReference("OutputTransform"))
    self.ui.outputSegmentationSelector.setCurrentNode(self._parameterNode.GetNodeReference("OutputSegmentationNode"))
    self.ui.checkBox_RotateAP.checked = (self._parameterNode.GetParameter("RotateAP") == "true")
    self.ui.inputVolumeSelector.setCurrentNode(self._parameterNode.GetNodeReference("InputVolume"))
    self.ui.outputVolumeSelector.setCurrentNode(self._parameterNode.GetNodeReference("OutputVolume"))

    # Update calculate button state and tooltip
    if self._parameterNode.GetNodeReference("InputSegmentationNode"):
      self.ui.calculateButton.toolTip = "Compute alignment angles"
      self.ui.calculateButton.enabled = True
    else:
      self.ui.calculateButton.toolTip = "Select a segmentation"
      self.ui.calculateButton.enabled = False

    # Show results if computed
    alignmentComputed = (self._parameterNode.GetParameter("OutputAlignmentAngleLR") != '' and \
                        (self._parameterNode.GetParameter("RotateAP") != "true" or self._parameterNode.GetParameter("OutputAlignmentAnglePA") != ''))
    if alignmentComputed:
      angleLR = float(self._parameterNode.GetParameter("OutputAlignmentAngleLR"))
      if self._parameterNode.GetParameter("RotateAP") == "true":
        anglePA = float(self._parameterNode.GetParameter("OutputAlignmentAnglePA"))
        paText = f', PA: {anglePA:.2f}°'
      else:
        paText = ''
      self.ui.resultsLabel.text = f"Alignment angles computed: LR: {angleLR:.2f}°{paText}"
    else:
      self.ui.resultsLabel.text = 'Alignment angles: Not yet calculated'

    # Enable/disable volume alignment section based on whether there are angles computed
    self.ui.alignVolumeCollapsibleButton.enabled = alignmentComputed

    # Update align volume button state and tooltip
    outputVolumeNode = self._parameterNode.GetNodeReference("OutputVolume")
    if alignmentComputed and self._parameterNode.GetNodeReference("InputVolume") and outputVolumeNode:
      self.ui.alignVolumeButton.toolTip = 'Crop and align input volume based on the computed alignment angles'
      self.ui.alignVolumeButton.enabled = True
    else:
      self.ui.alignVolumeButton.toolTip = 'Compute alignment angles and select input volume'
      self.ui.alignVolumeButton.enabled = False

    # Enable feeding input to analysis module if alignment is done for both the segment and the volume
    volumeAligned = False
    if outputVolumeNode:
      if outputVolumeNode.GetImageData() and outputVolumeNode.GetImageData().GetDimensions()[0] > 0:
        volumeAligned = True
    self.ui.feedInputToStrengthModuleButton.enabled = (alignmentComputed and volumeAligned)

    # All the GUI updates are done
    self._updatingGUIFromParameterNode = False

  def onInputChanged(self):
    # Reset results label text
    self.ui.resultsLabel.text = 'Alignment angles: Not yet calculated'

  def updateParameterNodeFromGUI(self, caller=None, event=None):
    """
    This method is called when the user makes any change in the GUI.
    The changes are saved into the parameter node (so that they are restored when the scene is saved and loaded).
    """

    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return

    wasModified = self._parameterNode.StartModify()  # Modify all properties in a single batch

    self._parameterNode.SetNodeReferenceID("InputSegmentationNode", self.ui.inputSegmentSelectorWidget.currentNodeID())
    self._parameterNode.SetParameter("InputSegmentID", str(self.ui.inputSegmentSelectorWidget.currentSegmentID()))
    self._parameterNode.SetNodeReferenceID("OutputTransform", self.ui.outputTransformSelectorComboBox.currentNodeID)
    self._parameterNode.SetNodeReferenceID("OutputSegmentationNode", self.ui.outputSegmentationSelector.currentNodeID)
    self._parameterNode.SetParameter("RotateAP", "true" if self.ui.checkBox_RotateAP.checked else "false")
    self._parameterNode.SetNodeReferenceID("InputVolume", self.ui.inputVolumeSelector.currentNodeID)
    self._parameterNode.SetNodeReferenceID("OutputVolume", self.ui.outputVolumeSelector.currentNodeID)

    self._parameterNode.EndModify(wasModified)

  def onCalculateButtonClicked(self):
    """
    Run processing when user clicks "Calculate" button.
    """
    # Check if output are empty, because otherwise it may result in data loss
    if (self._parameterNode.GetNodeReference("OutputTransform") and \
          self.logic.isTransformNodeNonIdentity(self._parameterNode.GetNodeReference("OutputTransform"))) or \
       (self._parameterNode.GetNodeReference("OutputSegmentationNode") and \
         self._parameterNode.GetNodeReference("OutputSegmentationNode").GetSegmentation().GetNumberOfSegments() > 0):
      result = slicer.util.confirmOkCancelDisplay('The output transform and/or the segmentation is non-empty!\n\n'
        'Any data in the transform will be overwritten, and the segmentation will be applied the new transform.\n\nDo you want to continue?',
        'Non-empty outputs! Data loss may occur')
      if not result:
        return  # User cancelled

    with slicer.util.tryWithErrorDisplay("Failed to compute results.", waitCursor=True):
      # Compute the alignment angles
      self.logic.calculateAlignmentAngles(self._parameterNode)
      self.updateGUIFromParameterNode()

  def onAlignVolumeButtonClicked(self):
    """
    Run processing when user clicks "Align volume" button.
    """
    with slicer.util.tryWithErrorDisplay("Failed to align the input volume.", waitCursor=True):
      # Crop and align the input volume using the calculated angles
      self.logic.cropAndAlignVolume(self._parameterNode)

  def onFeedInputToStrengthModuleButtonClicked(self):
    """
    Feed input to the Vertebral Strength Calculation module to facilitate analysis of aligned vertebrae.
    - Switch to the Vertebral Strength Calculation module
    - Set aligned volume as CT volume
    - Set segmentation containing the aligned segment
    - Create new strain and stress volumes
    """
    slicer.util.selectModule('VertebralStrengthCalculation')
    vscWidget = slicer.modules.vertebralstrengthcalculation.widgetRepresentation().self()

    alignedVolumeNode = self._parameterNode.GetNodeReference("OutputVolume")
    vscWidget.ui.inputCtSelector.setCurrentNode(alignedVolumeNode)

    alignedSegmentationNode = self._parameterNode.GetNodeReference("OutputSegmentationNode")
    vscWidget.ui.vertebraSegmentSelectorWidget.setCurrentNode(alignedSegmentationNode)
    segmentID = self._parameterNode.GetParameter("InputSegmentID")
    vscWidget.ui.vertebraSegmentSelectorWidget.setCurrentSegmentID(segmentID)

    vscWidget.ui.outputStressVolumeSelector.addNode()
    vscWidget.ui.outputStrainVolumeSelector.addNode()


#
# SegmentAxisAlignmentLogic
#

class SegmentAxisAlignmentLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self):
    """
    Called when the logic class is instantiated. Can be used for initializing member variables.
    """
    ScriptedLoadableModuleLogic.__init__(self)

    # set logic instance to the global variable that supplies it to the minimizer function
    global SegmentAxisAlignmentLogicInstanceGlobal
    SegmentAxisAlignmentLogicInstanceGlobal = self

    # Constants
    self.alignedVolumeNamePrefix = 'AlignedVolume_'

  def setDefaultParameters(self, parameterNode):
    """
    Initialize parameter node with default settings.
    """
    if not parameterNode.GetParameter("InitialAlignmentAngleLR"):
      parameterNode.SetParameter("InitialAlignmentAngleLR", '0')
    if not parameterNode.GetParameter("InitialAlignmentAnglePA"):
      parameterNode.SetParameter("InitialAlignmentAnglePA", '0')
    if not parameterNode.GetParameter("OutputAlignmentAngleLR"):
      parameterNode.SetParameter("OutputAlignmentAngleLR", '')
    if not parameterNode.GetParameter("OutputAlignmentAnglePA"):
      parameterNode.SetParameter("OutputAlignmentAnglePA", '')
    if not parameterNode.GetParameter("SegmentCroppingPaddingMm"):
      parameterNode.SetParameter("SegmentCroppingPaddingMm", '1')
    if not parameterNode.GetParameter("RotateAP"):
      parameterNode.SetParameter("RotateAP", 'true')

  def calculateAlignmentAngles(self, parameterNode):
    """
    Calculate angles for alignment corresponding to the minimum bounding box of the input segment.
    The resulting two angles (along LR and PA axes) are set as parameter into the parameter node.
    :param parameterNode: Parameter node containing the input
    """
    if not parameterNode:
      raise ValueError("Invalid parameter node given")
    inputSegmentationNode = parameterNode.GetNodeReference("InputSegmentationNode")
    inputSegmentID = parameterNode.GetParameter("InputSegmentID")
    if not inputSegmentationNode:
      raise ValueError("Invalid input segmentation node given")
    if inputSegmentID in [None, '']:
      raise ValueError("Invalid segment ID given")
    outputTransformNode = parameterNode.GetNodeReference("OutputTransform")
    if not outputTransformNode:
      raise ValueError("Invalid output transform node given")

    rotateAP = (parameterNode.GetParameter("RotateAP") == "true")

    import time
    startTime = time.time()
    logging.info('Finding alignment angles started')

    # Reset output parameters
    parameterNode.SetParameter("OutputAlignmentAngleLR", '')
    parameterNode.SetParameter("OutputAlignmentAnglePA", '')

    # Angles to be computed
    angleLR = 0.0
    anglePA = 0.0

    # Get closed surface representation poly data
    self.polygonSurface = vtk.vtkPolyData()
    inputSegmentationNode.CreateClosedSurfaceRepresentation()
    inputSegmentationNode.GetClosedSurfaceRepresentation(inputSegmentID, self.polygonSurface)

    #
    # Use optimizer to get the rotation angles corresponding to the
    # minimum volume bounding box of the convex hull
    #
    self.minimizer = vtk.vtkAmoebaMinimizer()
    self.minimizer.SetFunction(BoundingBoxMinimizerFunction)
    self.minimizer.SetParameterValue("rotLR", float(parameterNode.GetParameter("InitialAlignmentAngleLR")))
    self.minimizer.SetParameterScale("rotLR", 4)
    self.minimizer.SetParameterValue("rotPA", float(parameterNode.GetParameter("InitialAlignmentAnglePA")))
    if rotateAP:
      self.minimizer.SetParameterScale("rotPA", 4)
    else:
      self.minimizer.SetParameterScale("rotPA", 0)
    self.minimizer.SetMaxIterations(120)
    self.minimizer.Minimize()

    angleLR = self.minimizer.GetParameterValue("rotLR")
    if rotateAP:
      anglePA = self.minimizer.GetParameterValue("rotPA")

    self.setAlignmentTransformAngles(parameterNode, angleLR, anglePA)

    stopTime = time.time()
    logging.info(f'Finding alignment angles completed in {stopTime-startTime:.2f} seconds')

  def setAlignmentTransformAngles(self, parameterNode, angleLR, anglePA = 0.0):
    """
    Facilitate setting the alignment angles externally.
    """
    inputSegmentationNode = parameterNode.GetNodeReference("InputSegmentationNode")
    inputSegmentID = parameterNode.GetParameter("InputSegmentID")
    segmentName = inputSegmentationNode.GetSegmentation().GetSegment(inputSegmentID).GetName()
    if not inputSegmentationNode:
      raise ValueError("Invalid input segmentation node given")
    outputTransformNode = parameterNode.GetNodeReference("OutputTransform")
    if not outputTransformNode:
      raise ValueError("Invalid output transform node given")
    outputSegmentationNode = parameterNode.GetNodeReference("OutputSegmentationNode")

    alignmentTransform = vtk.vtkTransform()
    alignmentTransform.RotateZ(0.0)
    alignmentTransform.RotateY(anglePA)
    alignmentTransform.RotateX(angleLR)
    outputTransformNode.SetAndObserveTransformToParent(alignmentTransform)
    if self.isNodeNameDefault(outputTransformNode):
      outputTransformNode.SetName(f'AlignmentTransform_{inputSegmentationNode.GetName()}_{segmentName}')

    # Create new segmentation with selected segment and apply transform
    if outputSegmentationNode:
      if self.isNodeNameDefault(outputSegmentationNode):
        outputSegmentationNode.SetName(f'AlignedSegment_{inputSegmentationNode.GetName()}_{segmentName}')
      outputSegmentationNode.GetSegmentation().CopySegmentFromSegmentation(inputSegmentationNode.GetSegmentation(), inputSegmentID)
      outputSegmentationNode.CreateDefaultDisplayNodes()
      outputSegmentationNode.SetAndObserveTransformNodeID(outputTransformNode.GetID())

    # Set result in parameter node
    parameterNode.SetParameter("OutputAlignmentAngleLR", str(angleLR))
    parameterNode.SetParameter("OutputAlignmentAnglePA", str(anglePA))

  def cropAndAlignVolume(self, parameterNode):
    """
    - Crop input volume to the extent of the aligned segment with some padding
    - Supersample cropped volume to prevent data loss during rotation
    - Subsample the rotated volume (if requested, on by default)
    :param parameterNode: Parameter node containing the input
    """
    if not parameterNode:
      raise ValueError("Invalid parameter node given")
    outputSegmentationNode = parameterNode.GetNodeReference("OutputSegmentationNode")
    inputSegmentID = parameterNode.GetParameter("InputSegmentID")
    segmentName = outputSegmentationNode.GetSegmentation().GetSegment(inputSegmentID).GetName()
    if not outputSegmentationNode:
      raise ValueError("Invalid input segmentation node given")
    if inputSegmentID in [None, '']:
      raise ValueError("Invalid segment ID given")
    outputTransformNode = parameterNode.GetNodeReference("OutputTransform")
    if not outputTransformNode:
      raise ValueError("Invalid output transform node given")
    inputVolumeNode = parameterNode.GetNodeReference("InputVolume")
    if not inputVolumeNode:
      raise ValueError("Invalid input volume node given")
    outputVolumeNode = parameterNode.GetNodeReference("OutputVolume")
    if not outputVolumeNode:
      raise ValueError("Invalid output volume node given")

    # Disable the supersampling/subsampling steps, because the cropping module takes into account the rotation during cropping,
    # and as it is all done in one step, using interpolation for sampling the voxels to be rotated, no data loss occurs that we
    # could be able to mitigate with this supersampling/subsampling step.

    segmentCroppingPaddingMm = float(parameterNode.GetParameter("SegmentCroppingPaddingMm"))

    import time
    startTime = time.time()
    logging.info('Aligning volume started')

    #
    # Crop input volume to the extent of the input segment and supersample the cropped volume to prevent data loss during rotation
    #

    # Get aligned segment bounds
    outputSegmentationNode.CreateClosedSurfaceRepresentation()  # Make sure there is a poly data representation
    inputSegmentPolyData = outputSegmentationNode.GetClosedSurfaceInternalRepresentation(inputSegmentID)
    transformFilter = vtk.vtkTransformPolyDataFilter()
    transformFilter.SetTransform(outputTransformNode.GetTransformToParent())
    transformFilter.SetInputData(inputSegmentPolyData)
    transformFilter.Update()
    alignedSegmentPolyData = transformFilter.GetOutput()
    alignedSegmentBounds = [0]*6
    alignedSegmentPolyData.GetBounds(alignedSegmentBounds)

    # Add padding
    paddedBounds = [0]*6
    for i in range(3):
      paddedBounds[i*2] = alignedSegmentBounds[i*2] - segmentCroppingPaddingMm
      paddedBounds[i*2+1] = alignedSegmentBounds[i*2+1] + segmentCroppingPaddingMm

    # Create temporary ROI for cropping
    tempRoiNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsROINode')
    tempRoiNode.SetCenter((paddedBounds[0]+paddedBounds[1])/2, (paddedBounds[2]+paddedBounds[3])/2, (paddedBounds[4]+paddedBounds[5])/2)
    tempRoiNode.SetRadiusXYZ((paddedBounds[1]-paddedBounds[0])/2, (paddedBounds[3]-paddedBounds[2])/2, (paddedBounds[5]-paddedBounds[4])/2)

    # Get the highest level parent transform of the input volume. That is the one that needs to be transformed temporarily
    currentParentTransform = inputVolumeNode.GetParentTransformNode()
    nextParentTransform = currentParentTransform.GetParentTransformNode() if currentParentTransform else None
    while nextParentTransform is not None:
      currentParentTransform = nextParentTransform
      nextParentTransform = currentParentTransform.GetParentTransformNode()

    # Create cropped volume (so that the higher quality resample step can use it as a reference volume)
    tempReferenceVolumeNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLScalarVolumeNode', f'TempReferenceVolume')
    supersampleCropParams = slicer.vtkMRMLCropVolumeParametersNode()
    supersampleCropParams.SetInputVolumeNodeID(inputVolumeNode.GetID())
    supersampleCropParams.SetOutputVolumeNodeID(tempReferenceVolumeNode.GetID())
    supersampleCropParams.SetInterpolationMode(1)  # Use nearest neighbor to make it fast
    supersampleCropParams.SetROINodeID(tempRoiNode.GetID())
    slicer.mrmlScene.AddNode(supersampleCropParams)
    cropLogic = slicer.modules.cropvolume.logic()
    cropLogic.Apply(supersampleCropParams)
    slicer.mrmlScene.RemoveNode(supersampleCropParams)

    # Resample input volume using the transform
    resampleParameters = {'inputVolume': inputVolumeNode.GetID(), 'outputVolume': outputVolumeNode.GetID(),
      'referenceVolume': tempReferenceVolumeNode.GetID(), 'transformationFile': outputTransformNode.GetID(), 'interpolationType': 'bs'}
    slicer.cli.run(slicer.modules.resamplescalarvectordwivolume, None, resampleParameters, wait_for_completion=True)

    # Remove tempoarary ROI and reference volume
    slicer.mrmlScene.RemoveNode(tempRoiNode)
    slicer.mrmlScene.RemoveNode(tempReferenceVolumeNode)

    # Give meaningful name to output volume if it was the default
    if self.isNodeNameDefault(outputVolumeNode):
      outputVolumeNode.SetName(f'{self.alignedVolumeNamePrefix}{inputVolumeNode.GetName()}_{segmentName}')

    # Move transform and segmentation containing the aligned segment under the same parent as the aligned volume
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    outputVolumeItem = shNode.GetItemByDataNode(outputVolumeNode)
    parentItem = shNode.GetItemParent(outputVolumeItem)
    outputTransformItem = shNode.GetItemByDataNode(outputTransformNode)
    shNode.SetItemParent(outputTransformItem, parentItem)
    outputSegmentationItem = shNode.GetItemByDataNode(outputSegmentationNode)
    shNode.SetItemParent(outputSegmentationItem, parentItem)

    # Trigger UI update
    parameterNode.Modified()

    stopTime = time.time()
    logging.info(f'Aligning volume completed in {stopTime-startTime:.2f} seconds')

  def isNodeNameDefault(self, node):
    if not node:
      return True
    tag = node.GetNodeTagName()
    name = node.GetName()
    if name == tag:
      return True
    elif name[:len(tag)] == tag and name[len(tag)] == '_':
      try:
        int(name[len(tag)+1:])
        return True
      except ValueError:
        return False
    return False

  def isTransformNodeNonIdentity(self, transformNode):
    if not transformNode or not transformNode.GetTransformToParent():
      return False
    # Get the list of transforms in the composite transform
    isNonIdentityTransform = False
    transformsCollection = vtk.vtkCollection()
    slicer.vtkMRMLTransformNode.FlattenGeneralTransform(transformsCollection, transformNode.GetTransformToParent())
    identityMatrix = vtk.vtkMatrix4x4()
    for index in range(transformsCollection.GetNumberOfItems()):
      transform = transformsCollection.GetItemAsObject(index)
      # Linear: compare with slicer.vtkAddonMathUtilities.MatrixAreEqual() to an identity matrix
      if transform.IsA('vtkTransform'):
        isNonIdentityTransform |= not slicer.vtkAddonMathUtilities.MatrixAreEqual(transform.GetMatrix(), identityMatrix)
      # Grid and bspline: coefficient image values are all 0 (minimum and maximum scalar range is 0)
      elif transform.IsA('vtkGridTransform'):
        transform.Update()  # Make sure the displacement grid is available
        acc = vtk.vtkImageAccumulate()
        acc.SetInputData(transform.GetDisplacementGrid())
        acc.Update()
        isNonIdentityTransform |= max(acc.GetMax()) > 0
      else:
        logging.error('Unknown transform type %s' % (transform.GetClassName()))
    return isNonIdentityTransform

#
# Function to minimize the bounding box volume
#
def BoundingBoxMinimizerFunction():
  # get logic instance
  global SegmentAxisAlignmentLogicInstanceGlobal
  logic = SegmentAxisAlignmentLogicInstanceGlobal

  # get the current rotation parameters
  rotLR = logic.minimizer.GetParameterValue("rotLR")
  rotPA = logic.minimizer.GetParameterValue("rotPA")
  # logging.debug(f'Debugging: Iteration: {logic.minimizer.GetIterations()}: rotLR={rotLR:.2f}, rotPA={rotPA:.2f}')

  # compute the 4x4 transformation matrix for current angle values
  transform = vtk.vtkTransform()
  transform.RotateX(rotLR)
  transform.RotateY(rotPA)
  transform.RotateZ(0)
  transform.Update()

  transformFilter = vtk.vtkTransformPolyDataFilter()
  transformFilter.SetInputData(logic.polygonSurface)
  transformFilter.SetTransform(transform)

  # execute the 3D rotation of polygon surface
  transformFilter.Update()

  # get the bounding box volume of the rotated polygon surface
  transformedPolyData = transformFilter.GetOutput()
  bounds = transformedPolyData.GetBounds()
  volumeOfBoundingBox = (bounds[1] - bounds[0]) * (bounds[3] - bounds[2]) * (bounds[5] - bounds[4])

  logic.minimizer.SetFunctionValue(volumeOfBoundingBox)


#
# SegmentAxisAlignmentTest
#
class SegmentAxisAlignmentTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear()

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_SegmentAxisAlignment1()

  def test_SegmentAxisAlignment1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")

    self.delayDisplay('Loaded test data set')

    self.delayDisplay('Test passed')

