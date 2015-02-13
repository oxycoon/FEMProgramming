#include "window.h"

#include "glcontextsurfacewrapper.h"

//qt
#include <QQuickItem>
#include <QDebug>

Window::Window(QWindow *parent) :
  QQuickView(parent) {

  // Turn off clear before rendering !!!!!!!!!!!!!
  setClearBeforeRendering(false);

  // Set rendering surface type to GL surface
  setSurfaceType( OpenGLSurface );

  // Setup GLformat
  QSurfaceFormat fmt;
  fmt.setMajorVersion(4);
  fmt.setMinorVersion(3);
  fmt.setDepthBufferSize(24);
  fmt.setStencilBufferSize(1);
  fmt.setSwapBehavior( QSurfaceFormat::DoubleBuffer );
  fmt.setSamples(4);
  fmt.setProfile( QSurfaceFormat::CompatibilityProfile );
  setFormat(fmt);

  // Set the GLContext and scene graph to be persistent
  setPersistentOpenGLContext(true);
  setPersistentSceneGraph(true);

  // Initial size
  setResizeMode(QQuickView::SizeRootObjectToView);
  setMinimumSize( QSize( 800, 600 ) );

  // Create window
  create();
}

void Window::initGLSurface() {
  _glsurface = std::make_shared<GLContextSurfaceWrapper>(openglContext());
}
