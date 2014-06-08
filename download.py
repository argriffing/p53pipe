"""
Check this web page for new versions:
http://p53.fr/TP53_database_download/TP53_tumor_database/tumor_database.html

If the release has not changed,
The following text should appear verbatim in the html:
Current release: June 2012 (2012_R1)

Download the zipped text version of the disease database:
http://p53.fr/TP53_database_download/TP53_tumor_database/TP53_US/UMDTP53_curated_2012_R1_US.txt.zip

The extracted zip file should contain at least the following file:
UMDTP53_curated_2012_R1_US.txt

"""
from __future__ import division, print_function, absolute_import

from StringIO import StringIO
import urllib2
import zipfile


def main():

    html_url = '/'.join((
            'http://p53.fr',
            'TP53_database_download',
            'TP53_tumor_database',
            'tumor_database.html'))
    print('hardcoded website html url:', html_url)
    print('reading the html file from the internet...')
    response = urllib2.urlopen(html_url)
    html = response.read()
    print('successfully read the html file')
    print()

    current_release_text = 'Current release: June 2012 (2012_R1)'
    print('scanning the html text for the current release...')
    if current_release_text in html:
        print('found the expected current release')
    else:
        raise Exception('the site seems to have a new release')
    print()

    
    zip_url = '/'.join((
        'http://p53.fr',
        'TP53_database_download',
        'TP53_tumor_database',
        'TP53_US/UMDTP53_curated_2012_R1_US.txt.zip'))
    print('hardcoded data zip file url:', zip_url)
    print('reading the data zip file from the internet...')
    response = urllib2.urlopen(zip_url)
    z = response.read()
    print('successfully read the zip file')
    print()

    z_stringio = StringIO(z)
    data_filename = 'UMDTP53_curated_2012_R1_US.txt'
    zin = zipfile.ZipFile(z_stringio)
    print('extracting', data_filename, 'into the current directory...')
    zin.extract(data_filename)
    print('extracted the text data file')
    zin.close()
    print('closed the zip file')
    print()


main()

