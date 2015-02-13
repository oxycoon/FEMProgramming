#include "glscenerenderer.h"
#include "window.h"
#include "glcontextsurfacewrapper.h"

// Qt
#include <QQuickWindow>
#include <QDebug>

// stl
#include <iostream>




namespace Private {

  Renderer::Renderer (const std::string& name) : _name{name} {}

  Renderer::~Renderer() {}

  void
  Renderer::paint() {

//      qDebug() << "Painting on viewport " << _name.c_str() << ": " << _viewport;

    assert(_name.length());

//    qDebug() << "Renderer::paint; " << _name.c_str();

    const GMlib::TextureRenderTarget &render_tex = GMlibWrapper::getInstance().getRenderTextureOf(_name);
    const GMlib::GL::Texture& tex = render_tex.getTexture();

    if( !_prog.isValid() ) {
//      std::cout << "Prog ! valid: setting up." << std::endl;

      _vs.create();
      _fs.create();
      _prog.create();

      _vs.setSource(
            "#version 150 compatibility\n"
            "layout(std140) uniform;\n"
            "in  vec4 vertices;"
            "out vec2 coords;"
            "out vec4 gl_Position;"
            "void main() {"
            "    gl_Position = vertices;"
            "    coords = (vertices.xy + vec2(1.0,1.0))*0.5;"
            "}"
            );

      _fs.setSource(
            "#version 150 compatibility\n"
            "layout(std140) uniform;\n"
            "\n"
            "uniform sampler2D u_tex0;"
            "in vec2 coords;"
            "void main() {"
            "    gl_FragColor = texture( u_tex0, coords.st );"
            "}"
            );

      if( !_vs.compile() ) {
        std::cout << "Vertex shader compile error: " << _vs.getCompilerLog() << std::endl;
        exit(-666);
      }
      if(!_fs.compile()) {
        std::cout << "Fragment shader compile error: " << _fs.getCompilerLog() << std::endl;
        exit(-666);
      }
      _prog.attachShader(_vs);
      _prog.attachShader(_fs);
      if( _prog.link() != GL_TRUE ) {
        std::cout << "Render prog link error: " << _prog.getLinkerLog() << std::endl;
        exit(-666);
      }


      const GLfloat values[] = {
          -1.0f, -1.0f,
           1.0f, -1.0f,
          -1.0f,  1.0f,
           1.0f,  1.0f
      };

      _vbo.create();
      _vbo.bufferData( 2 * 4 * sizeof(GLfloat), values, GL_STATIC_DRAW );
    }


//      std::cout << "Prog OK!" << std::endl;


    glViewport(_viewport.x(), _viewport.y(), _viewport.width(), _viewport.height());

//      glDisable(GL_DEPTH_TEST);

//      glClearColor(0, 0, 0, 1);
//      glClear(GL_COLOR_BUFFER_BIT);

//      glDisable(GL_BLEND);
//      glEnable(GL_BLEND);
//      glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    _prog.bind(); {

//        _prog.setUniform( "u_tex0", _tex, GL_TEXTURE0, GLuint(0) );
      _prog.setUniform( "u_tex0", tex, GL_TEXTURE0, GLuint(0) );

      GMlib::GL::AttributeLocation vert_loc = _prog.getAttributeLocation("vertices");

      _vbo.bind();
      _vbo.enableVertexArrayPointer( vert_loc, 2, GL_FLOAT, GL_FALSE, 0, static_cast<const GLvoid*>(nullptr) );
      GL_CHECK(::glDrawArrays(GL_TRIANGLE_STRIP,0,4));
      _vbo.disable(vert_loc);
      _vbo.unbind();

    } _prog.unbind();


  }

  void
  Renderer::setName(const std::string &name) {
    _name = name;
  }

  void
  Renderer::setViewport(const QRectF &viewport ) {
    _viewport = viewport;
  }

  void
  Renderer::setViewport(int x, int y, int w, int h) {
    setViewport( QRectF(x,y,w,h) );
  }



}









GLSceneRenderer::GLSceneRenderer() : _renderer{nullptr}, _glsurface{nullptr}, _name{}, _paused{false} {

  _renderer.reset();
  setAcceptedMouseButtons(Qt::AllButtons);
  setFlag(ItemHasContents);
  setSmooth(false);
  connect( this, &QQuickItem::windowChanged, this, &GLSceneRenderer::handleWindowChanged );
}

GLSceneRenderer::~GLSceneRenderer() {

  _glsurface->makeCurrent(); {
    if(_renderer)
      _renderer.reset();
  } _glsurface->doneCurrent();
}

const
QString&  GLSceneRenderer::getTexName() const {

  return _name;
}

void
GLSceneRenderer::setTexName(const QString& name) {

  _name = name;

  if(!(_name.length() > 0))
    _renderer = nullptr;
//    _renderer.reset(nullptr);
  else if(_renderer)
    _renderer->setName(_name.toStdString());
}

bool
GLSceneRenderer::isPaused() const {

  return _paused;
}

void GLSceneRenderer::setPaused(bool paused) {

  _paused = paused;
}

void
GLSceneRenderer::sync() {

//  qDebug() << "_renderer" << _renderer.get();

  if(!_renderer && _name.length() > 0) {

    Window *w = static_cast<Window*>(window());
//    _renderer = std::unique_ptr<Private::Renderer>(new Private::Renderer(_name.toStdString()));
    _renderer = std::shared_ptr<Private::Renderer>(new Private::Renderer(_name.toStdString()));
    _glsurface = w->getGLSurface();
    connect( w, &Window::beforeRendering, _renderer.get(), &Private::Renderer::paint, Qt::DirectConnection );
  }

  if( !_renderer )
    return;

  QRectF r = mapRectToScene(boundingRect());
  QRectF cr(r.left(),window()->height()-r.bottom(), r.width(), r.height());
  _renderer->setViewport(cr);

  emit signViewportChanged(_name,r);
}

void
GLSceneRenderer::cleanup() {

  std::cout << "Cleanup" << std::endl;
}

void
GLSceneRenderer::handleWindowChanged(QQuickWindow* window) {

  if( !window ) return;

  Window *w = dynamic_cast<Window*>(window);
  if( !w ) return;

  connect( w, &Window::beforeSynchronizing, this, &GLSceneRenderer::sync );
  connect( w, &Window::sceneGraphInvalidated, this, &GLSceneRenderer::cleanup );
  connect( w, &Window::signFrameReady, this, &QQuickItem::update );
  connect( this, &GLSceneRenderer::signMouseMoved, w, &Window::signMouseMoved );
  connect( this, &GLSceneRenderer::signViewportChanged, w, &Window::signGuiViewportChanged );
  connect( this, &GLSceneRenderer::signMousePressed, w, &Window::signMousePressed );
  connect( this, &GLSceneRenderer::signMouseReleased, w, &Window::signMouseReleased);
  connect( this, &GLSceneRenderer::signMouseDoubleClicked, w, &Window::signMouseDoubleClicked);
  connect( this, &GLSceneRenderer::signKeyPressed, w, &Window::signKeyPressed );
  connect( this, &GLSceneRenderer::signKeyReleased, w, &Window::signKeyReleased);
  connect( this, &GLSceneRenderer::signWheelEventOccurred, w, &Window::signWheelEventOccurred);

}

void GLSceneRenderer::itemChange(ItemChange change, const ItemChangeData& value) {

//  qDebug() << "GLSceneRenderer changes: " << _name << ", change: " << change << ", value: " << value.boolValue;
  if(change == QQuickItem::ItemVisibleHasChanged && !value.boolValue)
    _renderer = nullptr;

  QQuickItem::itemChange(change,value);
}

QSGNode*GLSceneRenderer::updatePaintNode(QSGNode* node, QQuickItem::UpdatePaintNodeData*) {

//  qDebug() << "updatePaintNode; item has content: " << ((flags() & QQuickItem::ItemHasContents) == QQuickItem::ItemHasContents ? "yes" : "no");
  return node;
}

void GLSceneRenderer::mousePressEvent(QMouseEvent* event) {

  setFocus(true,Qt::MouseFocusReason);
  emit signMousePressed(_name,event);
}

void GLSceneRenderer::mouseReleaseEvent(QMouseEvent* event) {

  setFocus(true,Qt::MouseFocusReason);
  emit signMouseReleased(_name,event);
}

void GLSceneRenderer::keyPressEvent(QKeyEvent* event) {

  emit signKeyPressed(_name,event);
}

void GLSceneRenderer::keyReleaseEvent(QKeyEvent* event) {

  emit signKeyReleased(_name,event);
}

void GLSceneRenderer::wheelEvent(QWheelEvent* event) {

  setFocus(true,Qt::MouseFocusReason);
  emit signWheelEventOccurred(_name,event);
}

void GLSceneRenderer::mouseDoubleClickEvent(QMouseEvent* event) {

  setFocus(true,Qt::MouseFocusReason);
  emit signMouseDoubleClicked(_name,event);
}

void GLSceneRenderer::mouseMoveEvent(QMouseEvent *event)
{
    setFocus(true,Qt::MouseFocusReason);
    emit signMouseMoved(_name, event);
}

