#!/usr/bin/python3
#-*- encoding: utf-8 -*-

'''
read line database
'''


def rrl():
    """
    radio recombination lines

    Parameter
    =========
    None

    Return
    ======
    return: {mol:freq}
    """
    catalog = {}

    catalog["HI"] = [1420.40580]

    catalog["Halpha"] = [
    1013.76700, 1030.25100, 1047.09400, 1064.30700, 
    1081.89800, 1099.88000, 1118.26200, 1137.05600, 
    1156.27400, 1175.92700, 1196.02800, 1216.59000, 
    1237.62600, 1259.15000, 1281.17500, 1303.71800,
    1326.79200, 1350.41400, 1374.60000, 1399.36800, 
    1424.73400, 1450.71600, 1477.33500]
    
    catalog["Cbeta"] = [
    1001.14000, 1013.97500, 1027.03100, 1040.31200,
    1053.82200, 1067.56800, 1081.55400, 1095.78500,
    1110.26700, 1125.00500, 1140.00500, 1155.27300,
    1170.81500, 1186.63700, 1202.74600, 1219.14700,
    1235.84800, 1252.85500, 1270.17600, 1287.81800,
    1305.78800, 1324.09300, 1342.74300, 1361.74400,
    1381.10600, 1400.83600, 1420.94400, 1441.43900,
    1462.33000, 1483.62600
    ]
    
    catalog["Hebeta"] = [
    1001.04800, 1013.88200, 1026.93700, 1040.21600,
    1053.72600, 1067.47000, 1081.45500, 1095.68500,
    1110.16500, 1124.90200, 1139.90100, 1155.16800,
    1170.70800, 1186.52900, 1202.63600, 1219.03600,
    1235.73500, 1252.74100, 1270.06000, 1287.70000,
    1305.66800, 1323.97200, 1342.62000, 1361.62000,
    1380.97900, 1400.70800, 1420.81400, 1441.30700,
    1462.19600, 1483.49000
    ]

    return catalog


def from_cdms(db_file, Aij=None, Elow=None):
    """
    Get the line catalog from CDMS.

    Parameter
    ========
    db_file: the query output file from CDMS database.
    Aij: None | float. Set the lower Aij value.
    Elow: None | float. Set the upper Elow value.

    Return
    ======
    return: dict{mol:freq}
    """
    catalog = {}

    filter_out_mols = [
    "H2NCH2COOHI", "H2NCH2COOHII", "HCOCH2OH", "C2H3CN",
    "CH3OOH", "D2CO", "H2CO",
    "HCOOH", "cis-DCOOH", "DCOOH", "cis-HCOOD", "HCOOD"
    ]

    with open(db_file, "r") as fp:
        for line in fp.readlines():
            line = line.strip().split('|')

            if line[0] == "T_Name":
                continue

            mol = line[0].split(";")[0]
            rest_freq = float(line[3])

            if mol == "HCCCCCCCN":
                mol = "HC7N"
            if mol in filter_out_mols:
                continue

            if (Aij is not None) and (Aij < float(line[4])):
                continue

            if (Elow is not None) and (Elow < float(line[5])):
                continue

            catalog[mol] = [rest_freq]

    return catalog


def from_jpl(db_file):
    """
    Get the line catalog from JPL.

    Parameter
    ========
    db_file: the query output file from JPL database.

    Return
    ======
    return: dict{mol:freq}
    """
    catalog = {}

    filter_out_mols = ["C3H8O2"]

    with open(db_file, "r") as fp:
        for line in fp.readlines():
            line = line.strip().split()

            if len(line) < 5:
                mol = line[1]

                if line[0] == "75002":
                    mol = "H2NCH2COOHI"
                if line[0] == "75003":
                    mol = "H2NCH2COOHII"
                continue
            
            if mol in filter_out_mols:
                continue
            
            rest_freq = float(line[0][0:9])

            catalog[mol] = [rest_freq]

    return catalog


def from_splatalogue(db_file, Aij=None, Eup=None):
    """
    Get the line catalog from splatalogue query results.

    Parameter
    ========
    db_file: the query output file from splatalogue.
    Aij: None | float. Set the lower Aij value.
    Eup: None | float. Set the upper Eup value.

    Return
    ======
    return: dict{mol:freq}
    """
    catalog = {}

    with open(db_file, "r") as fp:
        for line in fp.readlines():
            line = line.strip().split(':')

            if line[0] == "Species":
                continue

            mol = line[0]
            rest_freq = float(line[2].split(',')[-1])
            
            if (Aij is not None) and (Aij < float(line[14])):
                continue

            if (Eup is not None) and (Eup < float(line[10])):
                continue

            catalog[mol] = [rest_freq]

    return catalog

