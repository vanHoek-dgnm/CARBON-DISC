import zipfile

zipfilePath = ("qcmon_1941-1950.zip")
zip = zipfile.ZipFile(zipfilePath)
zip.extractall(".")
zip.close()