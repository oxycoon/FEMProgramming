import QtQuick 2.2
import QtQuick.Controls 1.1
import QtQuick.Controls.Private 1.0 as Private
import QtQuick.Window 2.1


Item {
    id: root

    property int orientation: Qt.Horizontal
    property Component handleDelegate: Rectangle {
        width: 1
        height: 1
        color: Qt.darker(pal.window, 1.5)
    }

    property bool resizing: false

    /* \internal */
    default property alias __contents: contents.data
    /* \internal */
    property alias __items: splitterItems.children
    /* \internal */
    property alias __handles: splitterHandles.children

    clip: true
    Component.onCompleted: d.init()
    onWidthChanged: d.updateLayout()
    onHeightChanged: d.updateLayout()
    onOrientationChanged: d.changeOrientation()

    SystemPalette { id: pal }

    QtObject {
        id: d
        property bool horizontal: orientation == Qt.Horizontal
        readonly property string offset: horizontal ? "x" : "y"
        readonly property string otherOffset: horizontal ? "y" : "x"
        readonly property string size: horizontal ? "width" : "height"
        readonly property string otherSize: horizontal ? "height" : "width"
        readonly property string implicitSize: horizontal ? "implicitWidth" : "implicitHeight"
        readonly property string implicitOtherSize: horizontal ? "implicitHeight" : "implicitWidth"

        property int fillIndex: -1
        property bool updateLayoutGuard: true

        property bool validHandles : false


        function init()
        {

//          console.debug("d.init()")
//          console.debug("  root.size: " + root[d.size])
//          console.debug("  root.implicitSize: " + root[d.implicitSize])
//          console.debug("  items+handles size: " + splitterItems[d.size] + splitterHandles[d.size])
          var i
            for (i=0; i<__contents.length; ++i) {
                var item = __contents[i];
                if (!item.hasOwnProperty("x"))
                    continue

                if (splitterItems.children.length > 0)
                    handleLoader.createObject(splitterHandles, {"__handleIndex":splitterItems.children.length - 1})
                item.parent = splitterItems
                i-- // item was removed from list

                // should match disconnections in Component.onDestruction
                item.widthChanged.connect(d.updateLayout)
                item.heightChanged.connect(d.updateLayout)
            }
//          console.debug("  No. items: " + __items.length)
//          console.debug("  No. handles: " + __handles.length)

            for( i = 0; i < __handles.length; ++i)
              __handles[i].relativeOffset = i === 0 ? 1.0 / __items.length : __handles[i-1].relativeOffset + 1.0 / __items.length


            d.updateLayoutGuard = false
            d.updateLayout()
        }




        function changeOrientation()
        {
            if (__items.length == 0)
                return;
            d.updateLayoutGuard = true

            // Swap width/height for items and handles:
            for (var i=0; i<__items.length; ++i) {
                var item = __items[i]
                var tmp = item.x
                item.x = item.y
                item.y = tmp
                tmp = item.width
                item.width = item.height
                item.height = tmp

                var handle = __handles[i]
                if (handle) {
                    tmp = handle.x
                    handle.x = handle.y
                    handle.y = handle.x
                    tmp = handle.width
                    handle.width = handle.height
                    handle.height = tmp
                }
            }

            // Change d.horizontal explicit, since the binding will change too late:
            d.horizontal = orientation == Qt.Horizontal
            d.updateLayoutGuard = false
        }


        function updateLayout()
        {

//          console.debug("UPDATE LAYOUTS")

            // This function will reposition both handles and
            // items according to the their width/height:
            if (__items.length === 0)
                return;
            if (d.updateLayoutGuard === true)
                return
            d.updateLayoutGuard = true

//          console.debug("UPDATE LAYOUTS - PAST GUARDS")
//          console.debug("No. Handles: " + __handles.length )
//          console.debug("No. Items: " + __items.length )
          for( var i2 = 0; i2 < __handles.length; ++i2 )
//            console.debug("  handle[" + i2 + "] rel offset: " + __handles[i2].relativeOffset)

            // Position items and handles according to the relative handle positions
            var lastVisibleItem, lastVisibleHandle, handle, item
            for (var i=0; i<__items.length; ++i) {

                // Position item to the right of the previous visible handle:
                handle = __handles[i]
                if (handle && handle.visible) {
                    handle[d.offset] = handle.relativeOffset * root[d.size]
//                    console.debug("HANDLE OFFSET[" + i + "]: " + handle[d.offset])
                    handle[d.otherOffset] = 0
                    handle[d.otherSize] = root[d.otherSize]
                }

                item = __items[i];
                if( item.visible || i == d.fillIndex) {
                  item[d.offset] = lastVisibleHandle ? lastVisibleHandle[d.offset] + lastVisibleHandle[d.size] : 0
                  item[d.otherOffset] = 0
                  item[d.otherSize] = splitterItems[d.otherSize]

                  if(handle && handle.visible)
                    item[d.size] = lastVisibleHandle ? handle[d.offset] - lastVisibleHandle[d.offset] - lastVisibleHandle[d.size] : handle[d.offset]
                  else
                    item[d.size] = root[d.size] - lastVisibleHandle[d.offset] - lastVisibleHandle[d.size]


//                  console.debug("item [" + i + "]")
//                  console.debug("  s:  " + item[d.size])
//                  console.debug("  os: " + item[d.otherSize])


                  lastVisibleItem = item
                }
                lastVisibleHandle = handle

            }

            d.updateLayoutGuard = false
        }
    }

    Component {
        id: handleLoader
        Loader {
            id: itemHandle

            property int __handleIndex: -1
            property QtObject styleData: QtObject {
                readonly property int index: __handleIndex
                readonly property alias hovered: mouseArea.containsMouse
                readonly property alias pressed: mouseArea.pressed
                readonly property bool resizing: mouseArea.drag.active
                onResizingChanged: root.resizing = resizing
            }
            property real relativeOffset: 0.33
            property bool resizeLeftItem: (d.fillIndex > __handleIndex)
            visible: __items[__handleIndex].visible // __items[__handleIndex + (resizeLeftItem ? 0 : 1)].visible
            sourceComponent: handleDelegate
            onWidthChanged: d.updateLayout()
            onHeightChanged: d.updateLayout()
            onXChanged: { moveHandle(); d.updateLayout() }
            onYChanged: { moveHandle(); d.updateLayout() }

            MouseArea {
                id: mouseArea
                anchors.fill: parent
                property real defaultMargin: Private.Settings.hasTouchScreen ? Screen.pixelDensity * 3.5 : 2
                anchors.leftMargin: (parent.width <= 1) ? -defaultMargin : 0
                anchors.rightMargin: (parent.width <= 1) ? -defaultMargin : 0
                anchors.topMargin: (parent.height <= 1) ? -defaultMargin : 0
                anchors.bottomMargin: (parent.height <= 1) ? -defaultMargin : 0
                hoverEnabled: true
                drag.threshold: 0
                drag.target: parent
                drag.axis: root.orientation === Qt.Horizontal ? Drag.XAxis : Drag.YAxis
                cursorShape: root.orientation === Qt.Horizontal ? Qt.SplitHCursor : Qt.SplitVCursor
            }

            function moveHandle() {
                // Moving the handle means resizing an item. Which one,
                // left or right, depends on where the fillItem is.
                // 'updateLayout' will be overridden in case new width violates max/min.
                // 'updateLayout' will be triggered when an item changes width.
                if (d.updateLayoutGuard)
                    return

//                console.debug("Parent size: " + root[d.size])
//                console.debug(" handle offset: " + this[d.offset])
//                console.debug(" relative: " + this[d.offset]/root[d.size])
//                console.debug(" index: " + __handleIndex)

                relativeOffset = this[d.offset]/ Math.max(root[d.size],root[d.implicitSize])
            }

        }
    }

    Item {
        id: contents
        visible: false
        anchors.fill: parent
    }
    Item {
        id: splitterItems
        anchors.fill: parent
    }
    Item {
        id: splitterHandles
        anchors.fill: parent
    }

    Component.onDestruction: {
        for (var i=0; i<splitterItems.children.length; ++i) {
            var item = splitterItems.children[i];

            // should match connections in init()
            item.widthChanged.disconnect(d.updateLayout)
            item.heightChanged.disconnect(d.updateLayout)
        }
    }
}
