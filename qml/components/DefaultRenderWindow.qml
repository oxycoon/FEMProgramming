import QtQuick 2.1
import QtQuick.Layouts 1.1
import QtQuick.Controls 1.1

import "qrc:/qml/components"

Item {
  id: root

  View {
    id: singleview; anchors.fill: parent;
    visible: true;
    name: "Projection";
//    name: "Iso"
  }

//  RelativeSplitView {
//    id: multiview
//    anchors.fill: parent
//    orientation: Qt.Vertical

//    RelativeSplitView {
//      id: sp1
//      orientation: Qt.Horizontal
//      View {
//        id: top_view
//        visible: true; name: "Top";       onMouseDoubleClicked: root.state = "top_view";
//        Component.onCompleted: {
//          resizing = Qt.binding(d.viewsChanging)
//          keyPressed.connect(root.keyPressed)
//          keyReleased.connect(root.keyReleased)
//          mousePressed.connect(root.mousePressed)
//          mouseReleased.connect(root.mouseReleased)
//          mouseDoubleClicked.connect(root.mouseDoubleClicked)
//          mousePositionChanged.connect(root.mousePositionChanged)
//        }
//      }
//      View {
//        id: side_view
//        visible: true; name: "Side";      onMouseDoubleClicked: root.state = "side_view";
//        Component.onCompleted: {
//          resizing = Qt.binding(d.viewsChanging)
//          keyPressed.connect(root.keyPressed)
//          keyReleased.connect(root.keyReleased)
//          mousePressed.connect(root.mousePressed)
//          mouseReleased.connect(root.mouseReleased)
//          mouseDoubleClicked.connect(root.mouseDoubleClicked)
//          mousePositionChanged.connect(root.mousePositionChanged)
//        }
//      }
//    }

//    RelativeSplitView {
//      id: sp2
//      orientation: Qt.Horizontal
//      View {
//        id: projection_view
//        visible: true; name: "Projection"; onMouseDoubleClicked: root.state = "projection_view";
//        Component.onCompleted: {
//          resizing = Qt.binding(d.viewsChanging)
//          keyPressed.connect(root.keyPressed)
//          keyReleased.connect(root.keyReleased)
//          mousePressed.connect(root.mousePressed)
//          mouseReleased.connect(root.mouseReleased)
//          mouseDoubleClicked.connect(root.mouseDoubleClicked)
//          mousePositionChanged.connect(root.mousePositionChanged)
//        }
//      }
//      View {
//        id: front_view
//        visible: true; name: "Front";      onMouseDoubleClicked: root.state = "front_view";
//        Component.onCompleted: {
//          resizing = Qt.binding(d.viewsChanging)
//          keyPressed.connect(root.keyPressed)
//          keyReleased.connect(root.keyReleased)
//          mousePressed.connect(root.mousePressed)
//          mouseReleased.connect(root.mouseReleased)
//          mouseDoubleClicked.connect(root.mouseDoubleClicked)
//          mousePositionChanged.connect(root.mousePositionChanged)
//        }
//      }
//    }
//  }

  QtObject {
    id: d

    function viewsChanging() { return multiview.resizing || sp1.resizing || sp2.resizing }
  }

//  states: [
//    State{
//      name: "projection_view"
//      PropertyChanges { target: multiview; visible: false }
//      PropertyChanges { target: singleview; visible: true; name: projection_view.name }
//    },
//    State{
//      name: "top_view"
//      PropertyChanges { target: multiview; visible: false }
//      PropertyChanges { target: singleview; visible: true; name: top_view.name }
//    },
//    State{
//      name: "side_view"
//      PropertyChanges { target: multiview; visible: false }
//      PropertyChanges { target: singleview; visible: true; name: side_view.name }
//    },
//    State{
//      name: "front_view"
//      PropertyChanges { target: multiview; visible: false }
//      PropertyChanges { target: singleview; visible: true; name: front_view.name }
//    }
//  ]

}
