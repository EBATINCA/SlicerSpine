cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(SlicerSpine)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/SlicerSpine")
set(EXTENSION_CATEGORY "Spine")
set(EXTENSION_CONTRIBUTORS "Ron Alkalay (Beth Israel Deaconess Medical Center), Csaba Pinter (EBATINCA), Steve Pieper (Isomics)")
set(EXTENSION_DESCRIPTION "Tools useful in spine related research.")
set(EXTENSION_ICONURL "https://www.example.com/Slicer/Extensions/SlicerSpine.png")
set(EXTENSION_SCREENSHOTURLS "https://www.example.com/Slicer/Extensions/SlicerSpine/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
