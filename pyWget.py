import sys, urllib2

def saveDataFromURL(url):
    fileName = url.split('/')[-1]
    urlID    = urllib2.urlopen(url)
    fileID   = open(fileName,'wb')
    meta     = urlID.info()
    fileSize = int(meta.getheaders('Content-Length')[0])

    print 'Downloading: {0:s} Bytes: {1:d}'.format(fileName, fileSize)

    reaoutSize = 0
    blockSize  = 8192
    bar        = ''
    outString  = '\rBytes: {0:d} [{1:3.1f}%] Progress:[{2:>' + str(fileSize / blockSize) + 's}]'

    while True:
        bar += '='
        status = outString.format(reaoutSize, reaoutSize * 100. / fileSize, bar)
        sys.stdout.write(status)
        sys.stdout.flush()

        buffer = urlID.read(blockSize)
        if not buffer:
            break
        fileID.write(buffer)
        reaoutSize += len(buffer)

    print ''
    fileID.close()
