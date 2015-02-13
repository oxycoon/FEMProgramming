#ifndef __GLCONTEXTWRAPPER_H__
#define __GLCONTEXTWRAPPER_H__


// stl
#include <memory>


class QOpenGLContext;
class QOffscreenSurface;


class GLContextSurfaceWrapper {
public:
  explicit GLContextSurfaceWrapper(QOpenGLContext* context );

  void      makeCurrent();
  void      doneCurrent();

private:
  std::shared_ptr<QOpenGLContext>       _context;
  std::shared_ptr<QOffscreenSurface>    _surface;
};

#endif // __GLCONTEXTWRAPPER_H__
