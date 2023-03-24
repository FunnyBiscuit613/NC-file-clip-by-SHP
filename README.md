@[TOC]

# 前言
HELLO 介于之前发的《NC文件不规则裁剪(利用shp文件裁剪)》一文尚有些操作难度，所以本次打算分享一个更方便的NC文件不规则裁剪的方法。
# 示例数据
示例数据当然是我们的老演员NCEP数据了，下载链接：[点我下崽](https://downloads.psl.noaa.gov//Datasets/ncep.reanalysis/Dailies/surface_gauss/air.2m.gauss.2022.nc)

![在这里插入图片描述](https://img-blog.csdnimg.cn/bf1aae76c25146d7bb87d6d6085cebc0.png)
# 逻辑介绍
本次代码的主要逻辑就是让NC文件生成TIF，然后使用SHP裁剪TIF，最后再将TIF转为NC文件！
逻辑有点绕，不过好用就行！！！
# 特别说明
在上一版本的代码中会出现很多的维度问题，比如在试用具有level层的代码。大家可以先选出来自己需要的level，然后在进行。

# 代码
## 加载裤子

```python
from osgeo import gdal, ogr, osr
import numpy as np
import xarray as xr
import os
import shutil
```
## 打开nc文件
```python
# 打开nc文件
ds = xr.open_dataset(r'air.2m.gauss.2022.nc') # 打开nc文件
shp_path = r'D:/province.shp' # 设置SHP的路径
out_netCDF= r'air.2m.gauss.2022_clip.nc' # 生成的nc文件路径

lat = ds['lat'].values # 获取经纬度
lon = ds['lon'].values # 获取经纬度
air = ds['air'].values # 获取数据
time = ds['time'].values # 获取时间

```
## 生成网格信息

```python
# 生成网格信息
# 没有特殊需求的情况下，不需要修改这一段代码
LonMin,LatMax,LonMax,LatMin = [lon.min(),lat.max(),lon.max(),lat.min()]  # 获取经纬度范围
N_Lat = len(lat) # 获取经纬度个数
N_Lon = len(lon) # 获取经纬度个数
Lon_Res = (LonMax - LonMin) /(float(N_Lon)-1) # 获取经纬度分辨率
Lat_Res = (LatMax - LatMin) / (float(N_Lat)-1) # 获取经纬度分辨率
```

## 生成TIF

```python
# 创建一个临时文件夹
# 没有特殊需求的情况下，不需要修改这一段代码
'''
我也不知道为什呢，生成多波段的情况下数据会出错，所以我只能生成多个单波段的tif文件
如果有大佬知道怎么解决，欢迎留言
'''
if not os.path.exists('temporary_folders'): # 这里创建的临时文件夹也不用管，代码结束的时候会被删除
    os.makedirs('temporary_folders')

for i in range(len(air)):
    out_tif = 'temporary_folders/air.2m.gauss.2022_%s.tif' % i # 生成临时tif文件
    tif_ds = gdal.GetDriverByName('Gtiff').Create(out_tif,N_Lon,N_Lat,1,gdal.GDT_Float32) 
    geotransform = (LonMin,Lon_Res, 0, LatMin, 0, Lat_Res) # 设置栅格数据的左上角的经度坐标，经度分辨率，旋转角度，左上角的纬度坐标，纬度分辨率，旋转角度
    tif_ds.SetGeoTransform(geotransform) # 设置栅格数据的地理参考信息
    srs = osr.SpatialReference() # 创建坐标系
    srs.ImportFromEPSG(4326) # 定义地理坐标系统为WGS84(4326代表WGS84)
    tif_ds.SetProjection(srs.ExportToWkt()) # 设置栅格数据的坐标系
    tif_ds.GetRasterBand(1).WriteArray(air[i,:,:]) # 写入数据
    tif_ds.FlushCache() # 将数据写入文件
```
## 裁剪TIF

```python
temporary_folders_files = os.listdir('temporary_folders') # 获取临时文件夹下的文件
temporary_folders_files = [os.path.join('temporary_folders',i) for i in temporary_folders_files] # 获取临时文件夹下的文件路径


if not os.path.exists('temporary_folders_clip'):
    os.makedirs('temporary_folders_clip')

for i in range(len(temporary_folders_files)):
    ds = gdal.Open(temporary_folders_files[i]) # 打开刚刚输出的tif文件
    print(temporary_folders_files[i])
    print(r'temporary_folders_clip/air.2m.gauss.2022_clip_%s.tif' % i)
    
    # 使用shp裁剪
    ds = gdal.Warp(
    r'temporary_folders_clip/air.2m.gauss.2022_clip_%s.tif' % i, #待裁剪的影像完整路径（包括文件名）
    temporary_folders_files[i],    #裁剪后图像保存的完整路径（包括文件名）
    format='GTiff', # 保存图像的格式
    cutlineDSName=shp_path, # 矢量文件的完整路径
    cropToCutline=True, # 保证裁剪后影像大小跟矢量文件的图框大小一致（设置为False时，结果图像大小会跟待裁剪影像大小一样，则会出现大量的空值区域）
    # cutlineWhere=cutlineWhere, #矢量文件筛选条件,
    dstNodata=np.nan, #设置空值
    )
ds = None
```

##  重新生成NC文件

```python
# 重新生成NC文件

file_list = os.listdir('temporary_folders_clip') # 获取裁剪后的tif文件
file_list = [os.path.join('temporary_folders_clip',i) for i in file_list] # 获取裁剪后的tif文件路径
file_list

ds0 = gdal.Open(file_list[0]) # 打开第一个tif文件
data = ds0.ReadAsArray() # 读取数据
geotransform = ds0.GetGeoTransform() # 获取栅格数据的地理参考信息
num_x = ds0.RasterXSize # 获取栅格数据的列数
num_y = ds0.RasterYSize # 获取栅格数据的行数

print(num_x,num_y,'\n',geotransform)

# 生成网格信息
# 没有特殊需求的情况下，不需要修改这一段代码

lon = np.arange(geotransform[0],geotransform[0]+num_x*geotransform[1],geotransform[1]) # 生成经度

lat = np.arange(geotransform[3],geotransform[3]+num_y*geotransform[5],geotransform[5]) # 生成纬度

air_clip_nc = np.zeros((len(file_list),num_y,num_x)) # 生成一个空数组

# 获取数据
ds0_data = ds0.ReadAsArray() # 读取数据

# 获取array的序列
index_ = int(file_list[0].split('_')[-1].split('.')[0])
air_clip_nc[index_,:,:] = ds0_data

for i in range(1,len(file_list)):
    ds = gdal.Open(file_list[i]) # 打开tif文件
    ds_data = ds.ReadAsArray() # 读取数据
    index_ = int(file_list[i].split('_')[-1].split('.')[0])
    air_clip_nc[index_,:,:] = ds_data

# 生成nc文件
result_nc = xr.Dataset()
result_nc['air'] = xr.DataArray(air_clip_nc,coords=[('time',time),('lat',lat),('lon',lon)],dims=['time','lat','lon'])
result_nc['air'].attrs['long_name'] = 'air'
result_nc['air'].attrs['units'] = 'K'
result_nc['air'].attrs['standard_name'] = '2m_air_temperature'
result_nc['lon'].attrs['long_name'] = 'longitude'
result_nc['lon'].attrs['units'] = 'degrees_east'
result_nc['lon'].attrs['standard_name'] = 'longitude'
result_nc['lat'].attrs['long_name'] = 'latitude'
result_nc['lat'].attrs['units'] = 'degrees_north'
result_nc['lat'].attrs['standard_name'] = 'latitude'
result_nc['time'].attrs['long_name'] = 'time'

# 删除临时文件夹

shutil.rmtree('temporary_folders')
shutil.rmtree('temporary_folders_clip')

result_nc.to_netcdf(out_netCDF) # 保存nc文件

```

# 结果展示

![在这里插入图片描述](https://img-blog.csdnimg.cn/b9843cfd6d114c5fac13b18ed9f3a1f1.png)



看着效果还是很不错滴！！！

# 完整代码奉上

```python
from osgeo import gdal, ogr, osr
import numpy as np
import xarray as xr
import os
import shutil


ds = xr.open_dataset(r'air.2m.gauss.2022.nc') # 打开nc文件
out_netCDF = r'air.2m.gauss.2022_clip.nc' # 输出文件名
shp_path = r'D:\pythonGis/province.shp' # 获取shp文件的路径

lat = ds['lat'].values # 获取经纬度
lon = ds['lon'].values # 获取经纬度
air = ds['air'].values # 获取数据
time = ds['time'].values # 获取时间

# 生成网格信息
# 没有特殊需求的情况下，不需要修改这一段代码
LonMin,LatMax,LonMax,LatMin = [lon.min(),lat.max(),lon.max(),lat.min()]  # 获取经纬度范围
N_Lat = len(lat) # 获取经纬度个数
N_Lon = len(lon) # 获取经纬度个数
Lon_Res = (LonMax - LonMin) /(float(N_Lon)-1) # 获取经纬度分辨率
Lat_Res = (LatMax - LatMin) / (float(N_Lat)-1) # 获取经纬度分辨率

# 创建一个临时文件夹
# 没有特殊需求的情况下，不需要修改这一段代码
'''
我也不知道为什呢，生成多波段的情况下数据会出错，所以我只能生成多个单波段的tif文件
如果有大佬知道怎么解决，欢迎留言
'''
if not os.path.exists('temporary_folders'):
    os.makedirs('temporary_folders')

for i in range(len(air)):
    out_tif = 'temporary_folders/air.2m.gauss.2022_%s.tif' % i # 生成临时tif文件
    tif_ds = gdal.GetDriverByName('Gtiff').Create(out_tif,N_Lon,N_Lat,1,gdal.GDT_Float32) 
    geotransform = (LonMin,Lon_Res, 0, LatMin, 0, Lat_Res) # 设置栅格数据的左上角的经度坐标，经度分辨率，旋转角度，左上角的纬度坐标，纬度分辨率，旋转角度
    tif_ds.SetGeoTransform(geotransform) # 设置栅格数据的地理参考信息
    srs = osr.SpatialReference() # 创建坐标系
    srs.ImportFromEPSG(4326) # 定义地理坐标系统为WGS84(4326代表WGS84)
    tif_ds.SetProjection(srs.ExportToWkt()) # 设置栅格数据的坐标系
    tif_ds.GetRasterBand(1).WriteArray(air[i,:,:]) # 写入数据
    tif_ds.FlushCache() # 将数据写入文件
    
# 没有特殊需求的情况下，不需要修改这一段代码
temporary_folders_files = os.listdir('temporary_folders') # 获取临时文件夹下的文件
temporary_folders_files = [os.path.join('temporary_folders',i) for i in temporary_folders_files] # 获取临时文件夹下的文件路径

if not os.path.exists('temporary_folders_clip'):
    os.makedirs('temporary_folders_clip')

for i in range(len(temporary_folders_files)):
    ds = gdal.Open(temporary_folders_files[i]) # 打开刚刚输出的tif文件
    print(temporary_folders_files[i])
    print(r'temporary_folders_clip/air.2m.gauss.2022_clip_%s.tif' % i)
    
    # 使用shp裁剪
    ds = gdal.Warp(
    r'temporary_folders_clip/air.2m.gauss.2022_clip_%s.tif' % i, #待裁剪的影像完整路径（包括文件名）
    temporary_folders_files[i],    #裁剪后图像保存的完整路径（包括文件名）
    format='GTiff', # 保存图像的格式
    cutlineDSName=shp_path, # 矢量文件的完整路径
    cropToCutline=True, # 保证裁剪后影像大小跟矢量文件的图框大小一致（设置为False时，结果图像大小会跟待裁剪影像大小一样，则会出现大量的空值区域）
    # cutlineWhere=cutlineWhere, #矢量文件筛选条件,
    dstNodata=np.nan, #设置空值
    )
ds = None # 关闭文件
    
# 重新生成NC文件
file_list = os.listdir('temporary_folders_clip') # 获取裁剪后的tif文件
file_list = [os.path.join('temporary_folders_clip',i) for i in file_list] # 获取裁剪后的tif文件路径
file_list

ds0 = gdal.Open(file_list[0]) # 打开第一个tif文件
data = ds0.ReadAsArray() # 读取数据
geotransform = ds0.GetGeoTransform() # 获取栅格数据的地理参考信息
num_x = ds0.RasterXSize # 获取栅格数据的列数
num_y = ds0.RasterYSize # 获取栅格数据的行数

print(num_x,num_y,'\n',geotransform)

# 生成网格信息
# 没有特殊需求的情况下，不需要修改这一段代码

lon = np.arange(geotransform[0],geotransform[0]+num_x*geotransform[1],geotransform[1]) # 生成经度

lat = np.arange(geotransform[3],geotransform[3]+num_y*geotransform[5],geotransform[5]) # 生成纬度

air_clip_nc = np.zeros((len(file_list),num_y,num_x)) # 生成一个空数组

# 获取数据
ds0_data = ds0.ReadAsArray() # 读取数据

# 获取array的序列
index_ = int(file_list[0].split('_')[-1].split('.')[0])
air_clip_nc[index_,:,:] = ds0_data
for i in range(1,len(file_list)):
    ds = gdal.Open(file_list[i]) # 打开tif文件
    ds_data = ds.ReadAsArray() # 读取数据
    index_ = int(file_list[i].split('_')[-1].split('.')[0])
    air_clip_nc[index_,:,:] = ds_data

# 生成nc文件
result_nc = xr.Dataset()
result_nc['air'] = xr.DataArray(air_clip_nc,coords=[('time',time),('lat',lat),('lon',lon)],dims=['time','lat','lon'])
result_nc['air'].attrs['long_name'] = 'air'
result_nc['air'].attrs['units'] = 'K'
result_nc['air'].attrs['standard_name'] = '2m_air_temperature'
result_nc['lon'].attrs['long_name'] = 'longitude'
result_nc['lon'].attrs['units'] = 'degrees_east'
result_nc['lon'].attrs['standard_name'] = 'longitude'
result_nc['lat'].attrs['long_name'] = 'latitude'
result_nc['lat'].attrs['units'] = 'degrees_north'
result_nc['lat'].attrs['standard_name'] = 'latitude'
result_nc['time'].attrs['long_name'] = 'time'

# 删除临时文件夹

shutil.rmtree('temporary_folders')
shutil.rmtree('temporary_folders_clip')

result_nc.to_netcdf(out_netCDF) # 保存nc文件
```

 #  写在后面的话
上面的代码仅仅能处理三维的数据也就是(time, lat, lon)
如果您的数据是(level, time, lat, lon)这样的四维数据，大家可以先手动选择一下level
具体选择的方式如下：

```python
ds = xr.open_dataset(r'air.2m.gauss.2022.nc') # 打开nc文件
ds = ds.sel(level=1000) # 选择气压层
# 然后根据需要将数据的维度重新reshape成三维数据
# 加油！
```
