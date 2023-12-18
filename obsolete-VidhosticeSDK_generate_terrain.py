import os
import sys
sys.dont_write_bytecode = True		# disable creation of __pycache__
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))	# parent_dir
import VidhosticeSDK as sdk
os.system('title ' + __file__)

import re
import math
import requests
import json
from osgeo import gdal
# Stop GDAL printing both warnings and errors to STDERR
gdal.PushErrorHandler('CPLQuietErrorHandler')


name       = 'Vidhostice'
location   = '50.156N, 13.373E'
size       = 8192							# 8x8 km
angle      = 0								# pokud '' tak bude podle sjtsk jinak číslo (ve směru hod. ručiček)
resolution = 2048							# rozlišení heightmapy [maximálně 4096] (terrain.heightmap.png 1025 2049 4097)


def main():
	center = [0, 0]				# WGS84 (stupně) [Y, X] (Latitude, Longitude)
	pattern = '[\d]*[.][\d]+'
	i = 0
	if re.search(pattern, location) is not None:
		for catch in re.finditer(pattern, location):
			center[i] = float(catch[0])
			i += 1
	print('Střed mapy:', center)

	angle_sjtsk = -round((24.833332378-center[1])/1.34, 4)
	print('Pootočení křovákova zobrazení oproti WGS84:', angle_sjtsk)

	global angle
	if angle == '':
		angle = angle_sjtsk
	print('Otáčím mapu okolo středu o:', angle, 'stupňů (po směru hod. ručiček)\n')

	coords = []
	coords.append(center)
	# coordsX = [ vlevo nahoře, vpravo nahoře, vpravo dole, vlevo dole ]
	coords1 = getCoordsSet(center[0], center[1], size, angle)								# oblast co bude demgenem stažena
	coords2 = [ [max(coords1[0][0], coords1[1][0]), min(coords1[0][1], coords1[3][1])],		# bounding box všeho
				[max(coords1[0][0], coords1[1][0]), max(coords1[1][1], coords1[2][1])],
				[min(coords1[2][0], coords1[3][0]), max(coords1[1][1], coords1[2][1])],
				[min(coords1[2][0], coords1[3][0]), min(coords1[0][1], coords1[3][1])], ]
	coords3 = getCoordsSet(center[0], center[1], size // 2, angle)							# poloviční výřez (pouze pro převody 4km mapy na 2km)
	coords4 = [ [max(coords3[0][0], coords3[1][0]), min(coords3[0][1], coords3[3][1])],		# bounding box polovičního výřezu
				[max(coords3[0][0], coords3[1][0]), max(coords3[1][1], coords3[2][1])],
				[min(coords3[2][0], coords3[3][0]), max(coords3[1][1], coords3[2][1])],
				[min(coords3[2][0], coords3[3][0]), min(coords3[0][1], coords3[3][1])], ]
	coords = coords + coords1 + coords2 + coords3

	exportGPX(coords, filename='.mapy_cz.gpx', debug=False)

	with open('.downloadTiles.cmd', 'w', encoding='utf-8') as f:
		_mapycz_oblast = str(coords2[0][0])+' '+str(coords2[0][1])+' '+str(coords2[2][0])+' '+str(coords2[2][1])
		print('@echo off\nREM npm --location=global install sharp\nmkdir tiles\nREM export NODE_PATH=$(npm root -g)\nset NODE_PATH=C:\\Users\\teclv\\AppData\\Roaming\\npm\\node_modules\n\nset ZOOM=16\n', file=f)
		print('node mapycz-tile-downl.js tiles/ %ZOOM% bing         ' + _mapycz_oblast, file=f)
		print('REM                                                  Y vlevo nahore    X vlevo nahore     Y vpravo dole     X vpravo dole\n', file=f)
		print('REM bing base-m turist-m ophoto1618-m ophoto1415-m ophoto1012-m ophoto0406-m ophoto0203-m winter-m zemepis-m', file=f)
		print('REM army2-m (set ZOOM=15)\n', file=f)
		print('echo.\necho souradnice pro blender osm (blosm)', file=f)
		print('echo '+str(coords2[3][1])+','+str(coords2[3][0])+','+str(coords2[1][1])+','+str(coords2[1][0]), file=f)
		print('REM  X vlevo dole       Y vlevo dole       X vpravo nahore   Y vpravo nahore', file=f)
		print('echo.\n\necho.', file=f)
		print('REM https://www.openstreetmap.org/api/0.6/map?bbox='+str(coords4[3][1])+'%2C'+str(coords4[3][0])+'%2C'+str(coords4[1][1])+'%2C'+str(coords4[1][0]), file=f)
		print('echo.\n\necho.', file=f)
		print('echo Stred mapy:', center, '# WGS84 (stupne) [Y, X] (Latitude, Longitude)', file=f)
		print('echo Pootoceni krovakova zobrazeni oproti WGS84:', angle_sjtsk, file=f)
		print('echo Otacim mapu okolo stredu o:', angle, 'stupnu (po smeru hod. rucicek)', file=f)
		print('echo.\n\npause', file=f)


	demgenDir = os.path.dirname(os.path.realpath(__file__)).replace('\\', '/') + '/demgenDir'
	try: 
		os.mkdir(demgenDir)
		print('Vytvářím složku:', demgenDir)
	except OSError as error:
		print('Složka již existuje:', demgenDir)

	with open('.demgen.xml', 'w', encoding='utf-8') as f:
		print('<?xml version=\'1.0\' encoding=\'utf-8\'?>\n<demgenConfig>\n  <country value="CZ"/>', file=f)
		print('  <centralPoint latitude="' + str(center[0]) + '" longitude="' + str(center[1]) + '"/>', file=f)
		print('  <terrain size="' + str(size) + '" metersPerPixel="' + str(size / resolution) + '" rotation="' + str(angle) + '"/>', file=f)
		print('  <dataSavingFolder path="' + demgenDir + '"/>', file=f)
		print('</demgenConfig>', file=f)


#	exportCUZK(coords, download=False, link=True, filename='.cuzk_cz')
#	exportCUZK(coords, download=True, link=False, filename='.cuzk_cz')

#	# zde nastavit center buď na nulu nebo ručně podle potřeby
#	print('\n### zde nastavit center buď na nulu nebo ručně podle potřeby ###')
#	center = 127
#	print('Center:', center)
#	minHeight, maxHeight, heightScale = gdal_info('.cuzk_cz.tiff', center, verbose=False)
#	print()

#	gdal_convert('.cuzk_cz.tiff', '.terrain.heightmap_' + str(resolution) + 'x' + str(resolution) + '_heightScale_' + str(heightScale) + '.png', minHeight, maxHeight)


### ------------------------------------------------------------------------------------------------------ ###


def exportGPX(coords, filename='.mapy_cz.gpx', debug=False):
	with open(filename, 'w', encoding='utf-8') as f:
		print('<?xml version="1.0" encoding="utf-8"?><gpx xmlns="http://www.topografix.com/GPX/1/1" version="1.1" creator="VidhosticeSDK">', file=f)

		if not debug:
			for j in range(1, 30, 4):	# až 8 tracků/čtverců za sebou
				if len(coords) > j:
					print('<trk><name>' + str(j) + ' / okraj</name><trkseg>', file=f)
					for i in range(j, j+4):
						print('<trkpt lat="' + str(coords[i][0]) + '" lon="' + str(coords[i][1]) + '"></trkpt>', file=f)
					print('<trkpt lat="' + str(coords[j][0]) + '" lon="' + str(coords[j][1]) + '"></trkpt>', file=f)
					print('</trkseg></trk>', file=f)

		print('<wpt lat="' + str(coords[0][0]) + '" lon="' + str(coords[0][1]) + '"><name>0 / střed mapy</name></wpt>', file=f)
		if debug:
			for k in range(1, len(coords)):
				print('<wpt lat="' + str(coords[k][0]) + '" lon="' + str(coords[k][1]) + '"><name>' + str(k) + '</name></wpt>', file=f)

		print('</gpx>', file=f)
		print('Soubor "' + filename + '" vytvořen pro načtení v prohlížeči na stránce mapy.cz\n')


def gdal_info(inputImg, center=0, verbose=True):
	if verbose:
		sdk.info('Analyzuji soubor: ' + inputImg + '\n\n' + gdal.Info(inputImg, deserialize=True))
	else:
		sdk.info('Analyzuji soubor: ' + inputImg)

	gtif = gdal.Open(inputImg)
	srcband = gtif.GetRasterBand(1)
	stats = srcband.ComputeStatistics(0)
	minHeight = int(stats[0])
	maxHeight = int(stats[1])
	meanHeight = int(stats[2])
	rozdil = maxHeight - minHeight
	pulka = rozdil // 2
	print('Minimální nadmořská výška:', minHeight)
	print('Maximální nadmořská výška:', maxHeight)
	print('Průměrná nadmořská výška: ', meanHeight)
	print('---------------------------------')
	if rozdil < 255:
		print('Rozdíl max. - min. výšek: ', rozdil, '(méně jak heightScale="255" v map.i3d - OK)')
	else:
		sdk.error('Rozdíl max. - min. výšek:  ', rozdil, 'POZOR - bude potřeba poladit hodnoty')

	out_heightScale = 255

	if center == 0:
		stred = minHeight + pulka
		print('\nStřed vypočítán automaticky na:', stred)
		out_minHeight = stred - 127
		out_maxHeight = stred + 128
	else:
		stred = center
		print('\nStřed zadán ručně na:', stred)
		out_minHeight = stred - 127
		out_maxHeight = stred + 128

	print('Nastavuji MIN:', out_minHeight)
	print('Nastavuji MAX:', out_maxHeight)

	return out_minHeight, out_maxHeight, out_heightScale


def gdal_convert(inputImg, outputImg, minHeight, maxHeight):
	sdk.info('Konvertuji vstupní obrázek TIFF do 16ti bitového černobílého obrázku PNG')
	print('Ukládám "' + outputImg + '" - ', end='', flush=True)
	#
	# https://gis.stackexchange.com/questions/352643/gdal-translate-in-python-where-do-i-find-how-to-convert-the-command-line-argum
	# gdal_translate.exe -of PNG -ot UInt16 -scale %min% %max% 0 65535 ..exportImage.tif map_dem.png
	#
	kwargs = {
		'format': 'PNG',
		'outputType': gdal.GDT_UInt16,
		'scaleParams': [[minHeight, maxHeight, 0, 65535]]
	}
	ds = gdal.Translate(outputImg, inputImg, **kwargs)
	# do something with ds if you need
	ds = None # close and save ds
	print('Hotovo')


def exportCUZK(coords, download=False, link=False, filename='.cuzk_cz'):
	x1, y1 = wgs84_to_sjtsk_loc(coords[1][0], coords[1][1])
	x2, y2 = wgs84_to_sjtsk_loc(coords[3][0], coords[3][1])

	if download:
		sdk.info('Stahuji geodata z cuzk.cz')
		sdk.downloadURL('https://ags.cuzk.cz/arcgis2/rest/services/dmr5g/ImageServer/exportImage?bbox=' +
			str(x1) + '%2C' + str(y1) + '%2C' + str(x2) + '%2C' + str(y2) +
			'&bboxSR=&size=' + str(resolution) + '%2C' + str(resolution) + '&imageSR=&time=&format=tiff&pixelType=F32' +
			'&noData=&noDataInterpretation=esriNoDataMatchAny&interpolation=+RSP_BilinearInterpolation&compression=&compressionQuality=' +
			'&bandIds=&sliceId=&mosaicRule=&renderingRule=&adjustAspectRatio=true&validateExtent=false&lercVersion=1&compressionTolerance=&f=image', filename + '.tiff')
		sdk.downloadURL('https://ags.cuzk.cz/arcgis2/rest/services/dmr5g/ImageServer/exportImage?bbox=' +
			str(x1) + '%2C' + str(y1) + '%2C' + str(x2) + '%2C' + str(y2) +
			'&bboxSR=&size=' + str(resolution) + '%2C' + str(resolution) + '&imageSR=&time=&format=png&pixelType=F32' +
			'&noData=&noDataInterpretation=esriNoDataMatchAny&interpolation=+RSP_BilinearInterpolation&compression=&compressionQuality=' +
			'&bandIds=&sliceId=&mosaicRule=&renderingRule=&adjustAspectRatio=true&lercVersion=1&compressionTolerance=&f=image', filename + '.png')
	if link:
		sdk.info('Vytvářím linky pro geodata z cuzk.cz')
		with open(filename + '_TIFF.url', 'w', encoding='utf-8') as f:
			print('[InternetShortcut]', file=f)
			print('URL=https://ags.cuzk.cz/arcgis2/rest/services/dmr5g/ImageServer/exportImage?bbox=' +
			str(x1) + '%2C' + str(y1) + '%2C' + str(x2) + '%2C' + str(y2) +
			'&bboxSR=&size=' + str(resolution) + '%2C' + str(resolution) + '&imageSR=&time=&format=tiff&pixelType=F32' +
			'&noData=&noDataInterpretation=esriNoDataMatchAny&interpolation=+RSP_BilinearInterpolation&compression=&compressionQuality=' +
			'&bandIds=&sliceId=&mosaicRule=&renderingRule=&adjustAspectRatio=true&validateExtent=false&lercVersion=1&compressionTolerance=&f=image', file=f)
		with open(filename + '_PNG.url', 'w', encoding='utf-8') as f:
			print('[InternetShortcut]', file=f)
			print('URL=https://ags.cuzk.cz/arcgis2/rest/services/dmr5g/ImageServer/exportImage?bbox=' +
			str(x1) + '%2C' + str(y1) + '%2C' + str(x2) + '%2C' + str(y2) +
			'&bboxSR=&size=' + str(resolution) + '%2C' + str(resolution) + '&imageSR=&time=&format=png&pixelType=F32' +
			'&noData=&noDataInterpretation=esriNoDataMatchAny&interpolation=+RSP_BilinearInterpolation&compression=&compressionQuality=' +
			'&bandIds=&sliceId=&mosaicRule=&renderingRule=&adjustAspectRatio=true&lercVersion=1&compressionTolerance=&f=html', file=f)


def wgs84_to_sjtsk_net(B, L, H = 89.79):
	res = requests.get('http://epsg.io/trans?x='+str(L)+'&y='+str(B)+'&s_srs=4326&t_srs=5514&ops=1623')
	sjtsk_net = json.loads(res.text)

	for i in sjtsk_net:
		sjtsk_net[i] = float(sjtsk_net[i])

	return sjtsk_net['x'], sjtsk_net['y']


# Original JS code by Tomáš Pecina (tomas@pecina.cz), source: https://www.pecina.cz/krovak.html
# To Python3 translated by crpl – 2019-11-14
def wgs84_to_sjtsk_loc(B, L, H = 89.79):
	d2r = math.pi / 180
	a = 6378137.0
	f1 = 298.257223563
	dx = -570.69
	dy = -85.69
	dz = -462.84
	wx = 4.99821 / 3600 * math.pi / 180
	wy = 1.58676 / 3600 * math.pi / 180
	wz = 5.2611 / 3600 * math.pi / 180
	m  = -3.543e-6

	B *= d2r
	L *= d2r

	e2 = 1 - math.pow(1 - 1 / f1, 2)
	rho = a / math.sqrt(1 - e2 * math.pow(math.sin(B), 2))
	x1 = (rho + H) * math.cos(B) * math.cos(L)
	y1 = (rho + H) * math.cos(B) * math.sin(L)
	z1 = ((1 - e2) * rho + H) * math.sin(B)

	x2 = dx + (1 + m) * (x1 + wz * y1 - wy * z1)
	y2 = dy + (1 + m) * (-wz * x1 + y1 + wx * z1)
	z2 = dz + (1 + m) * (wy * x1 - wx * y1 + z1)

	a = 6377397.15508
	f1 = 299.152812853

	ab = f1 / (f1 - 1)
	p = math.sqrt(math.pow(x2, 2) + math.pow(y2, 2))
	e2 = 1 - math.pow(1 - 1 / f1, 2)
	th = math.atan(z2 * ab / p)
	st = math.sin(th)
	ct = math.cos(th)
	t = (z2 + e2 * ab * a * math.pow(st, 3)) / (p - e2 * a * math.pow(ct, 3))

	B = math.atan(t)
	H = math.sqrt(1 + math.pow(t, 2)) * (p - a / math.sqrt(1 + (1 - e2) * math.pow(t, 2)))
	L = 2 * math.atan(y2 / (p + x2))

	a = 6377397.15508
	e = 0.081696831215303
	n = 0.97992470462083
	rho0 = 12310230.12797036
	sinUQ = 0.863499969506341
	cosUQ = 0.504348889819882
	sinVQ = 0.420215144586493
	cosVQ = 0.907424504992097
	alpha  = 1.000597498371542
	k2 = 1.00685001861538

	sinB = math.sin(B)
	t = (1 - e * sinB) / (1 + e * sinB)
	t = math.pow(1 + sinB, 2) / (1 - math.pow(sinB, 2)) * math.exp(e * math.log(t))
	t = k2 * math.exp(alpha * math.log(t))

	sinU = (t - 1) / (t + 1)
	cosU = math.sqrt(1 - math.pow(sinU, 2))
	V = alpha * L
	sinV = math.sin(V)
	cosV = math.cos(V)
	cosDV = cosVQ * cosV + sinVQ * sinV
	sinDV = sinVQ * cosV - cosVQ * sinV
	sinS = sinUQ * sinU + cosUQ * cosU * cosDV
	cosS = math.sqrt(1 - math.pow(sinS, 2))
	sinD = sinDV * cosU / cosS
	cosD = math.sqrt(1 - math.pow(sinD, 2))

	eps = n * math.atan(sinD / cosD)
	rho = rho0 * math.exp(-n * math.log((1 + sinS) / cosS))

	CX = rho * math.sin(eps)
	CY = rho * math.cos(eps)
	return -CX, -CY


# map corners coordinates calculating
def getCoordsSet(b, l, dim, angle):
	brad = math.radians(b)
	lrad = math.radians(l)
	dist1 = dim * math.sqrt(2) * math.cos(math.radians(45 - math.fabs(angle))) / 2
	dist2 = dim * math.sqrt(2) * math.sin(math.radians(45 - math.fabs(angle))) / 2

	# https://en.wikipedia.org/wiki/Latitude#Length_of_a_degree_of_latitude
	# (phi - 0.5 deg) → (phi + 0.5 deg)
	latlen = 111132.954 - 559.822 * math.cos(2 * brad) + 1.175 * math.cos(4 * brad)
	# https://en.wikipedia.org/wiki/Longitude#Length_of_a_degree_of_longitude
	def lonlen(lat):
		return math.radians(math.pi / math.radians(180) * 6378137 * math.cos(math.atan(6356752.3 / 6378137 * math.tan(lat))))

	boffset1 = math.radians(dist1 / latlen)
	boffset2 = math.radians(dist2 / latlen)
	loffset1 = math.radians(dist1 / lonlen(brad)) # centralny punkt / central point
	loffset2 = math.radians(dist2 / lonlen(brad + boffset1)) # północ / north
	loffset3 = math.radians(dist2 / lonlen(brad - boffset1)) # południe / south

	# [NW, NE, SE, SW]
	if angle > 0:
		coordsSet = [[math.degrees(brad + boffset1), math.degrees(lrad - loffset2)], [math.degrees(brad + boffset2), math.degrees(lrad + loffset1)], [math.degrees(brad - boffset1), math.degrees(lrad + loffset3)], [math.degrees(brad - boffset2), math.degrees(lrad - loffset1)]]
	else:
		coordsSet = [[math.degrees(brad + boffset2), math.degrees(lrad - loffset1)], [math.degrees(brad + boffset1), math.degrees(lrad + loffset2)], [math.degrees(brad - boffset2), math.degrees(lrad + loffset1)], [math.degrees(brad - boffset1), math.degrees(lrad - loffset3)]]

	return coordsSet


if __name__ == '__main__':
	sdk.welcome()
	sdk.info('Generuji potřebné soubory pro vytvoření výškové mapy.')
	main()
	sdk.key_press('HOTOVO')
