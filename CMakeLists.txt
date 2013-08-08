cmake_minimum_required(VERSION 2.8.9)

project(TransformVisualizer)

# ------------------------------------------------------------------------
set(TRANSFORMVISUALIZER_VERSION_MAJOR "0")
set(TRANSFORMVISUALIZER_VERSION_MINOR "1")
set(TRANSFORMVISUALIZER_VERSION_PATCH "0")
set(TRANSFORMVISUALIZER_VERSION ${TRANSFORMVISUALIZER_VERSION_MAJOR}.${TRANSFORMVISUALIZER_VERSION_MINOR}.${TRANSFORMVISUALIZER_VERSION_PATCH})

#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/TransformVisualizer";)
set(EXTENSION_CATEGORY "Registration")
set(EXTENSION_CONTRIBUTORS "Franklin King (PerkLab, Queen's University), Csaba Pinter (PerkLab, Queen's University), Andras Lasso (PerkLab, Queen's University)")
set(EXTENSION_DESCRIPTION " Extension version: ${TRANSFORMVISUALIZER_VERSION}.")
set(EXTENSION_ICONURL "http://www.slicer.org/slicerWiki/images/5/55/TransformVisualizerLogo.png";)
set(EXTENSION_SCREENSHOTURLS "http://www.slicer.org/slicerWiki/images/f/f7/SlicerRT_DefFieldVisGlyphUI.png http://www.slicer.org/slicerWiki/images/a/aa/SlicerRT_DefFieldVisModels.png";)

#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
add_subdirectory(TransformVisualizer)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})