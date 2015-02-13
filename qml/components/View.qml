import QtQuick 2.2
import QtQuick.Controls 1.1

import MyCustomQuick 1.0


Item {
  id: root

  property string name : ""
  property bool   resizing : false

  signal forceRender

  Rectangle {
    id: bah
    color: "gray"
    anchors.fill: parent
    visible: false // root.resizing
  }

  GLSceneRenderer{
    id: renderer
    anchors.fill: parent
    name: root.name
//      visible: !root.resizing
    paused: false
//      view_type: root.view_type

    Label {
      anchors.top: parent.top
      anchors.left: parent.left
      anchors.margins: 5
      text: root.name
      color: "white"
      font {
        bold: true
        weight: Font.Black
      }
    }
  }

  onResizingChanged: {
    d.has_been_resized = d.has_been_resized | resizing

    if( !resizing && d.has_been_resized ) {
      d.has_been_resized = false
      root.forceRender()
    }
  }

  onVisibleChanged: forceRender()

  QtObject {
    id: d
    property bool   has_been_resized : false
  }

//  Component.onCompleted: root.forceRender.connect(renderer.forceRender)
}
