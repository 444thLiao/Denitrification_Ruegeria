"""
Docstring for figures.ForQGIS.partialCodeForQGIS
"""


name = '12_22Dry'
indir = "XXXXXXXXX"
output_tif = '{indir}/{name}DIN.tif'.format(name=name)
processing.run("qgis:idwinterpolation", {'INTERPOLATION_DATA':'file://{indir}/{name}.tsv?type=csv&delimiter=%5Ct&maxFields=10000&detectTypes=yes&xField=x&yField=y&crs=EPSG:4326&spatialIndex=no&subsetIndex=no&watchFile=no::~::0::~::32::~::0'.format(name=name),'DISTANCE_COEFFICIENT':2,'EXTENT':'113.745000000,114.438000000,22.000600000,22.590600000 [EPSG:4326]','PIXEL_SIZE':0.0005,'OUTPUT':output_tif})



# Following code are how to generate a suitable color bar.
## you need to use following code in the python console of QGIS 3.44
rlayer = iface.addRasterLayer(output_tif, name)
r = rlayer.renderer()
r.type()
provider = rlayer.dataProvider()
stats = provider.bandStatistics(1, QgsRasterBandStats.All)
if 'AllDIN.tif' not in output_tif:
    val_min=3.588088274002075
    val_max=204.99058532714844
    val_mid=104.28933680057526
else:
    val_min = stats.minimumValue
    val_max = stats.maximumValue
    val_mid = (val_min + val_max) / 2
print(val_mid)
# 定义颜色（根据你的图）
col_min = QColor(255, 245, 225)   # 浅黄
col_mid = QColor(255, 180, 120)   # 橙色
col_max = QColor(180, 30, 0)      # 深红

# 创建连续渐变色带
shader = QgsColorRampShader()
shader.setColorRampType(QgsColorRampShader.Interpolated)

items = [
    QgsColorRampShader.ColorRampItem(val_min, col_min, f"{val_min:.2f}"),
    QgsColorRampShader.ColorRampItem(val_mid, col_mid, f"{val_mid:.2f}"),
    QgsColorRampShader.ColorRampItem(val_max, col_max, f"{val_max:.2f}"),
]

shader.setColorRampItemList(items)

# 创建渲染器
fnc = QgsRasterShader()
fnc.setRasterShaderFunction(shader)

renderer = QgsSingleBandPseudoColorRenderer(rlayer.dataProvider(), 1, fnc)

# 强制图例显示 μM 范围，而不是 0~255
renderer.setClassificationMin(val_min)
renderer.setClassificationMax(val_max)

rlayer.setRenderer(renderer)
rlayer.triggerRepaint()
print(val_min,val_max,val_mid)
print(stats.minimumValue,stats.maximumValue,(stats.minimumValue + stats.maximumValue) / 2)

### after doing this, you need to manually transfer the color bar from continuous to quantile. This make the color represeting lower values more distinct.

# Following code are for making a print layout in qgis

# 获取当前工程
project = QgsProject.instance()

# 创建 Print Layout 管理器
manager = project.layoutManager()

# 创建一个 layout
layout = QgsPrintLayout(project)
layout.initializeDefaults()
layout.setName(name)

# 添加 layout 到工程
manager.addLayout(layout)

# -------------- 设置地图项 ----------------------

# 创建地图项
map_item = QgsLayoutItemMap(layout)

# 设置地图项的位置和大小（单位：毫米，可改）
map_item.setRect(20, 20, 200, 150)

# 指定你要显示的范围 extent
extent = QgsRectangle(113.745000000, 22.000600000,
                      114.438000000, 22.590600000)

# 设置地图项的显示范围
map_item.setExtent(extent)

# 把地图项加到 layout 里
layout.addLayoutItem(map_item)

# 强制刷新 layout 视图
map_item.attemptResize(QgsLayoutSize(200, 150))
map_item.refresh()