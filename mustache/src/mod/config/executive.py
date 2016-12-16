import sys, os
import argparseTypes
import json
import wget
import gzip
import tarfile
from pprint import pprint

def action(data_dir):
    config_file = os.path.join(data_dir,'mustache.config')

    config = {}

    if is_config(config_file):
        print "Config File Contents: "
        pprint(dict(json.load(open(config_file))))
        sys.exit()
    else:

        if download_taxondb():
            taxondb = os.path.abspath(os.path.join(data_dir, 'TaxonomyDatabase'))
            os.makedirs(taxondb)
            filename = wget.download("ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", taxondb)
            print
            print "Extracting the downloaded TaxonomyDatabase..."
            untar(filename, True)
            config['taxondb'] = taxondb

        else:

            taxonomyDB_loc = get_taxonomyDB_location()
            config['taxondb'] = taxonomyDB_loc


    print "Writing mustache.config file..."
    json.dump(config, open(config_file, 'w'), indent=4, sort_keys=True)
    print "Done :)"


def is_config(config_file):
    if os.path.isfile(config_file):
        print "Config file found at %s" % config_file
        return True
    else:
        print "No config file found, creating one at %s" % config_file
        return False


def download_taxondb():
    response = raw_input("An NCBI taxonomy database folder has not been specified. Would you like to download one? (y/n) --> ")
    if response.lower()[0] == 'y':
        return True
    else:
        return False


def get_taxonomyDB_location():
    taxonDB = ''
    response = raw_input("Then please give the path to your local NCBI taxonomy database --> ").strip()
    while True:
        try:
            taxonDB = argparseTypes.taxonDB(response)
            break
        except:
            response = raw_input("The path you gave does not contain a valid NCBI taxonomy database. Please try again --> ").strip()
            continue
    return os.path.abspath(taxonDB)


def untar(fname, delete=False):
    directory = os.path.dirname(fname)

    if (fname.endswith("tar.gz")):
        tar = tarfile.open(fname)
        tar.extractall(directory)
        tar.close()

        if delete:
            os.remove(fname)
    else:
        print "Not a tar.gz file: '%s '" % fname