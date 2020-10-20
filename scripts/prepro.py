import re
import os
import sys
import numpy as np
from pathlib import Path
#
MT = {18 : 'fission', 16 : '(n,2n)', 17 : '(n,3n)',
      102 : '(n,gamma)', 103 : '(n,p)', 107 : '(n,a)'}

def parse_prepros(path, numproc=1):
    """Parse prepro results

    Parametres:
    -----------
    path : str
        path to PREPRO calculation results
    numproc : int
        process number

    Returns:
    --------
    nuclidpath : iterable of str
        path to nuclides files
    nuclidnames : iterable of str
        name of nuclides"""
    nuclidpath = []
    nuclidnames = []
    for i in range(numproc):
        fg = [g for g in Path(os.path.join(path, str(i))).glob('*.tab')]
        nuclidpath = nuclidpath + fg
        for f in fg:
            nuclidnames.append(f.name.split(f.suffix))
    return nuclidpath, nuclidnames


def prepare_xs(path, numbergroup=1):
    """Prepare the needed representation of cross-section data

     Paramteres:
     -----------
     path : str
        filename of cross-section data
     numbergroup : int
        number of energies neutron multigroup

     Returns:
     --------
     energies : iterable of str
        energy discritization by multigroups
     xslib : dict
        key : MT number, value cross-section values (str)
     """
    def skip(ofile, number):
        for i in range(number):
            line = next(ofile)
    energies = np.zeros(numbergroup + 1)
    xslib = {}
    xs = []
    mtnum = ''
    with open(path,'r') as f:
        for line in f:
            res = re.search("MT=\w*\d+", line)
            if res:
               mtnum = re.search("\d+", line).group()
               skip(f, 5)
               xs = np.zeros(numbergroup)
               while(len(line.rstrip()) > 1):
                   dump = line.rstrip().split()
                   num = 0
                   en = 0.0
                   x = 0.0
                   for i, d in enumerate(dump):
                       if (i % 3 == 0):
                           num = int(d.rstrip())
                       if (i % 3 == 1):
                           en = float(d.rstrip())
                           if (num < numbergroup + 2):
                               if (energies[num - 1] == 0.0):
                                   energies[num - 1] = en
                       if (i % 3 == 2):
                           x = float(d.rstrip())
                           if (num < numbergroup + 1):
                               xs[num - 1] = x
                   line = next(f)
               if (sum(xs) > 0):
                   xslib[mtnum] = xs
    return energies, xslib

def create_xml(nuclidpath, nuclidnames, numbergroup, flux = None):
    """Creating a xsmg.xml file with cross-section data

     Paramteres:
     -----------
     nuclidpath : iterable of str
        filename of cross-section data
     nuclidnames : iterable of str
        nuclide names
     numbergroup : int
        number of energies neutron multigroup
     flux : iterable of double (optional)
        weight function to collapse cross-section into reaction rates
    """
    result = {}
    energies = None
    for p, n in zip(nuclidpath, nuclidnames):
        energies, result[n] = prepare_xs(p)
    from lxml import etree
    root = etree.Element("impxslib")
    childe = etree.SubElement(root, "energy", ng=str(numbergroup))
    childe.text = ' '.join(["{:10.4e}".format(e) for e in energies[:]])
    childl = etree.SubElement(root, "xslibs", typex="xs", ng=str(numbergroup))
    for k, v in result.items():
        for kk, vv in v.items():
            if (int(kk) in MT.keys()):
                childx = etree.SubElement(childl, "xslib", ng=str(numbergroup),
                 reaction = MT[int(kk)], name = k)
                #childx.text = "{:10.4e}".format((interflux * v.sf).sum())
                childx.text = ' '.join(["{:10.4e}".format(e) for e in vv])
    tree=etree.ElementTree(root)
    tree.write('xsmg.xml', xml_declaration=True, encoding='utf-8',
     pretty_print=True)

if __name__ == '__main__':
    numbergroup = int(sys.argv[-1])
    num = int(sys.argv[-2])
    gpath = sys.argv[-3]
    nuclidpath, nuclidnames = parse_prepros(gpath, num, numbergroup)
    create_xml(nuclidpath, nuclidnames, numbergroup)