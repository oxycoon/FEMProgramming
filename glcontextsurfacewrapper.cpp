#include "glcontextsurfacewrapper.h"

#include "utils.h"

// qt
#include <QOpenGLContext>
#include <QOffscreenSurface>
#include <QDebug>

// stl
#include <stdexcept>

GLContextSurfaceWrapper::GLContextSurfaceWrapper(QOpenGLContext* context) {

  if( !context->isValid() ) throw std::invalid_argument("OpenGLContext provided not valid!" + __EXCEPTION_TAIL);

  qDebug() << "Initialized GL Context Format: " << context->format();

  // Create Internal shared GL context
  _context = std::make_shared<QOpenGLContext>();
  _context->setShareContext( context );
  _context->setFormat( context->format() );
  _context->create();

  // Set up offscreen rendering surface for GMlib rendering
  _surface = std::make_shared<QOffscreenSurface>();
  _surface->setFormat(_context->format());
  _surface->create();
}

void GLContextSurfaceWrapper::makeCurrent() {

  _context->makeCurrent(_surface.get());
}

void GLContextSurfaceWrapper::doneCurrent() {

  _context->doneCurrent();
}

