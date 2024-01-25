import numpy as np
import pandas as pd

def loadEBSD_crc(fname, *args, **kwargs):
    def getCS(cpr):
        CS = []
        for p in range(1, cpr['phases']['count'] + 1):
            phase = cpr['phases']['phase' + str(p)]

            if 'spacegroup' in phase and phase['spacegroup'] > 0:
                Laue = ('spaceId', phase['spacegroup'])
            elif 'lauegroup' in phase and isinstance(phase['lauegroup'], str):
                Laue = [phase['lauegroup']]
            else:
                LaueGroups = ['-1', '2/m', 'mmm', '4/m', '4/mmm', '-3', '-3m', '6/m', '6/mmm', 'm-3', 'm-3m']
                Laue = [LaueGroups[phase['lauegroup']]]
            
            CS.append(crystalSymmetry(Laue, [phase['a'], phase['b'], phase['c']],
                                      [phase['alpha'], phase['beta'], phase['gamma']] * degree,
                                      'mineral', phase['structurename']))
        CS = ['notIndexed'] + CS
        return CS

    def getJobParam(cpr):
        job = cpr['job']
        param = {'cells': False, 'unitCell': None, 'm': None, 'n': None, 'x': None, 'y': None,
                 'ColumnNames': ['Phase'], 'ColumnType': [1]}

        ColumnNames = [
            'X', 'Y', 'phi1', 'Phi', 'phi2', 'MAD', 'BC', 'BS', 'Unknown', 'Bands', 'Error', 'ReliabilityIndex'
        ]
        dataType = [4] * 6 + [1] * 5 + [4]

        for k in range(1, cpr['fields']['count'] + 1):
            order = cpr['fields']['field' + str(k)]

            if order <= 12:
                param['ColumnNames'].append(ColumnNames[order - 1])
                param['ColumnType'].append(dataType[order - 1])
            else:
                param['ColumnNames'].append('Unknown' + str(order))
                param['ColumnType'].append(4)

        return param

    def localCRCLoader(crcFile, params):
        with open(crcFile, 'rb') as file:
            data = np.fromfile(file, dtype=np.uint8)

        if not params['n']:
            n = len(data) // params['m']
        else:
            n = params['n']

        data = data.reshape(-1, n)
        type = params['ColumnType']
        ndx = np.cumsum([0] + type)

        d = np.zeros((len(type), n), dtype=np.float32)
        d[type == 4, :] = data[ndx[type == 4][:, None] + np.arange(4)]
        d[type != 4, :] = data[1 + ndx[type != 4][:, None]]

        if params['cells']:
            d = np.vstack((d, params['x'], params['y'])
            params['ColumnNames'].extend(['x', 'y'])

        data_dict = {params['ColumnNames'][i]: d[i] for i in range(len(params['ColumnNames'])}
        loader = pd.DataFrame(data_dict)

        return loader

    def localCPRParser(cprFile):
        with open(cprFile, 'rb') as file:
            str = file.read().decode('utf-8')

        cpr = {}
        Title = None
        lines = str.strip().split('\n')
        for line in lines:
            line = line.strip()

            if line.startswith('['):
                Title = line[1:-1].replace(' ', '').lower()
                cpr[Title] = {}
            elif line:
                field, value = line.split('=')
                field = field.replace(' ', '').lower()
                value = value.strip()

                try:
                    numericValue = float(value)
                    cpr[Title][field] = numericValue
                except ValueError:
                    cpr[Title][field] = value

        return cpr

    try:
        if not fname.endswith(('.cpr', '.crc')):
            raise ValueError("File extension must be .cpr or .crc")

        path, file = os.path.split(fname)
        cprFile = os.path.join(path, file + '.cpr')
        crcFile = os.path.join(path, file + '.crc')

        cpr = localCPRParser(cprFile)

        CS = kwargs.get('CS', getCS(cpr))
        param = getJobParam(cpr)

        if 'check' in kwargs:
            ebsd = EBSD()
            return ebsd

        loader = localCRCLoader(crcFile, param)
        q = loader['getRotations']()
        phases = loader['getColumnData']('Phase')
        options = loader['getOptions']('ignoreColumns', 'phase')

        ebsd = EBSD(q, phases, CS, options, 'unitCell', param['unitCell'])
        ebsd.opt.cprInfo = cpr

    except Exception as e:
        interfaceError(fname)

    if 'convertSpatial2EulerReferenceFrame' in kwargs:
        ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, 180 * degree), 'keepEuler')
    elif 'convertEuler2SpatialReferenceFrame' in kwargs:
        ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, 180 * degree), 'keepXY')
    elif 'wizard' not in kwargs:
        warning(
            ".crc and .cpr files usually have inconsistent conventions for spatial coordinates and Euler angles. "
            "You may want to use one of the options 'convertSpatial2EulerReferenceFrame' or 'convertEuler2SpatialReferenceFrame' to correct for this")
    
    return ebsd
