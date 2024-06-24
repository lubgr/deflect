local bvp = import 'bvp.libsonnet';
local test = import 'test.libsonnet';

{
  // This boundary value problem is published on https://labs.degreetutors.com, which links against
  // the spreadsheet
  // https://docs.google.com/spreadsheets/d/1FmxJPeD_yA1AKOaWs_-vYf5IIV2LSzM9oZJAuyBxsNo/
  name: 'power_pole',

  material: bvp.LinElast('default', E=200000e6, nu=0.3, rho=78000 / 9.81),
  crosssection: bvp.Generic('default', A=0.05, Iyy=10e-6, Izz=10e-6),

  // Avoid typing out numeric keys all the time:
  local indexPlusOneToKey(array) = {
    [std.toString(i + 1)]: array[i]
    for i in std.range(0, std.length(array) - 1)
  },

  nodes: indexPlusOneToKey([
    [-4.0000, -4.0000, 0.0000],
    [-1.3361, 1.3361, 26.0473],
    [-1.3361, 1.3361, 21.1204],
    [4.0000, -4.0000, 0.0000],
    [1.3361, -1.3361, 19.0000],
    [-1.3361, 1.3361, 30.9118],
    [-1.9069, 1.9069, 14.9289],
    [3.2372, -3.2372, 5.4403],
    [-2.7145, -2.7145, 9.1685],
    [-2.7145, 2.7145, 9.1685],
    [2.7145, -2.7145, 9.1685],
    [-2.3259, -2.3259, 11.9405],
    [-2.1138, 2.1138, 13.4530],
    [2.3259, 2.3259, 11.9405],
    [2.3259, -2.3259, 11.9405],
    [2.1138, -2.1138, 13.4530],
    [1.3361, 1.3361, 28.6666],
    [-1.3361, 1.3361, 28.6666],
    [1.9069, 1.9069, 14.9289],
    [1.9069, -1.9069, 14.9289],
    [-1.7225, -1.7225, 16.2442],
    [1.7225, 1.7225, 16.2442],
    [1.7225, -1.7225, 16.2442],
    [-1.5490, -1.5490, 17.4813],
    [1.5490, -1.5490, 17.4813],
    [1.3361, 1.3361, 23.8645],
    [1.3361, -1.3361, 30.9118],
    [1.3361, 1.3361, 30.9118],
    [1.3361, -1.3361, 21.1204],
    [1.3361, 1.3361, 21.1204],
    [1.5490, 1.5490, 17.4813],
    [1.3361, -1.3361, 23.8645],
    [-1.5490, 1.5490, 17.4813],
    [1.3361, -1.3361, 26.0473],
    [-1.7225, 1.7225, 16.2442],
    [1.3361, 1.3361, 26.0473],
    [1.3361, -1.3361, 28.6666],
    [2.9156, -0.5528, 24.9258],
    [-1.3361, 1.3361, 19.0000],
    [-3.2372, 3.2372, 5.4403],
    [-3.2372, -3.2372, 5.4403],
    [2.1138, 2.1138, 13.4530],
    [-2.3259, 2.3259, 11.9405],
    [-2.1138, -2.1138, 13.4530],
    [-1.9069, -1.9069, 14.9289],
    [2.7145, 2.7145, 9.1685],
    [3.2372, 3.2372, 5.4403],
    [1.3361, 1.3361, 19.0000],
    [5.0526, -0.5528, 24.1205],
    [-1.3361, -1.3361, 26.0473],
    [4.0000, 4.0000, 0.0000],
    [-4.0000, 4.0000, 0.0000],
    [0.0000, -3.2372, 5.4403],
    [3.2372, 0.0000, 5.4403],
    [0.0000, 3.2372, 5.4403],
    [-3.2372, 0.0000, 5.4403],
    [0.0000, 0.0000, 32.3389],
    [2.9156, 0.5528, 20.0220],
    [5.0526, 0.5528, 19.2168],
    [5.0526, -0.5528, 19.2168],
    [2.9156, -0.5528, 20.0220],
    [-1.3361, -1.3361, 19.0000],
    [-1.3361, -1.3361, 21.1204],
    [5.0526, 0.5528, 24.1205],
    [2.9156, 0.5528, 24.9258],
    [-2.9156, 0.5528, 20.0220],
    [-5.0526, 0.5528, 19.2168],
    [-5.0526, -0.5528, 19.2168],
    [-2.9156, -0.5528, 20.0220],
    [2.9156, -0.5528, 29.6216],
    [5.0526, -0.5528, 28.8163],
    [5.0526, 0.5528, 28.8163],
    [2.9156, 0.5528, 29.6216],
    [-1.3361, 1.3361, 23.8645],
    [-1.3361, -1.3361, 30.9118],
    [-1.3361, -1.3361, 23.8645],
    [-1.3361, -1.3361, 28.6666],
    [-2.9156, -0.5528, 24.9258],
    [-5.0526, -0.5528, 24.1205],
    [-5.0526, 0.5528, 24.1205],
    [-2.9156, 0.5528, 24.9258],
    [-2.9156, -0.5528, 29.6216],
    [-5.0526, -0.5528, 28.8163],
    [-5.0526, 0.5528, 28.8163],
    [-2.9156, 0.5528, 29.6216],
  ]),

  elements: indexPlusOneToKey(std.map(function(nodes) bvp.Truss3d(nodes=nodes), [
    ['8', '4'],
    ['47', '55'],
    ['11', '8'],
    ['15', '11'],
    ['11', '9'],
    ['9', '10'],
    ['9', '12'],
    ['16', '15'],
    ['15', '12'],
    ['14', '15'],
    ['20', '16'],
    ['22', '19'],
    ['23', '20'],
    ['19', '20'],
    ['25', '23'],
    ['23', '21'],
    ['22', '23'],
    ['21', '24'],
    ['5', '25'],
    ['25', '24'],
    ['37', '27'],
    ['6', '17'],
    ['5', '29'],
    ['30', '29'],
    ['2', '17'],
    ['26', '2'],
    ['29', '32'],
    ['3', '26'],
    ['48', '3'],
    ['32', '34'],
    ['42', '7'],
    ['36', '34'],
    ['2', '18'],
    ['34', '37'],
    ['9', '15'],
    ['11', '12'],
    ['12', '16'],
    ['18', '6'],
    ['28', '6'],
    ['20', '21'],
    ['21', '25'],
    ['23', '24'],
    ['7', '35'],
    ['5', '24'],
    ['42', '13'],
    ['36', '37'],
    ['32', '36'],
    ['30', '32'],
    ['5', '30'],
    ['22', '25'],
    ['20', '22'],
    ['19', '23'],
    ['16', '19'],
    ['14', '16'],
    ['11', '14'],
    ['27', '70'],
    ['7', '13'],
    ['56', '41'],
    ['28', '18'],
    ['48', '33'],
    ['35', '31'],
    ['45', '44'],
    ['17', '18'],
    ['17', '72'],
    ['37', '71'],
    ['17', '28'],
    ['31', '33'],
    ['48', '31'],
    ['35', '33'],
    ['26', '65'],
    ['46', '47'],
    ['43', '44'],
    ['51', '54'],
    ['47', '51'],
    ['32', '38'],
    ['34', '65'],
    ['1', '41'],
    ['10', '40'],
    ['53', '41'],
    ['41', '9'],
    ['43', '10'],
    ['43', '14'],
    ['12', '43'],
    ['12', '44'],
    ['16', '44'],
    ['20', '45'],
    ['45', '21'],
    ['36', '38'],
    ['34', '50'],
    ['15', '44'],
    ['44', '20'],
    ['16', '45'],
    ['45', '23'],
    ['36', '65'],
    ['34', '38'],
    ['32', '64'],
    ['26', '49'],
    ['32', '50'],
    ['50', '37'],
    ['39', '24'],
    ['10', '12'],
    ['9', '43'],
    ['4', '54'],
    ['48', '5'],
    ['54', '8'],
    ['14', '46'],
    ['46', '11'],
    ['10', '46'],
    ['42', '14'],
    ['19', '42'],
    ['42', '16'],
    ['31', '22'],
    ['35', '22'],
    ['21', '35'],
    ['31', '25'],
    ['28', '27'],
    ['30', '48'],
    ['26', '30'],
    ['36', '26'],
    ['26', '32'],
    ['39', '30'],
    ['39', '31'],
    ['7', '22'],
    ['17', '36'],
    ['13', '19'],
    ['14', '13'],
    ['36', '2'],
    ['17', '37'],
    ['30', '3'],
    ['19', '7'],
    ['43', '13'],
    ['55', '40'],
    ['52', '40'],
    ['48', '39'],
    ['33', '39'],
    ['28', '37'],
    ['27', '17'],
    ['34', '17'],
    ['26', '34'],
    ['29', '26'],
    ['48', '29'],
    ['48', '25'],
    ['5', '31'],
    ['23', '31'],
    ['42', '20'],
    ['15', '42'],
    ['46', '15'],
    ['39', '3'],
    ['36', '18'],
    ['22', '33'],
    ['13', '44'],
    ['19', '35'],
    ['43', '42'],
    ['7', '45'],
    ['46', '43'],
    ['10', '14'],
    ['26', '64'],
    ['2', '50'],
    ['32', '49'],
    ['73', '72'],
    ['71', '70'],
    ['44', '7'],
    ['13', '45'],
    ['47', '54'],
    ['24', '33'],
    ['50', '18'],
    ['21', '33'],
    ['35', '24'],
    ['45', '35'],
    ['7', '21'],
    ['12', '13'],
    ['8', '53'],
    ['1', '53'],
    ['4', '53'],
    ['52', '55'],
    ['51', '55'],
    ['40', '56'],
    ['52', '56'],
    ['1', '56'],
    ['9', '53'],
    ['11', '53'],
    ['11', '54'],
    ['46', '54'],
    ['10', '55'],
    ['46', '55'],
    ['9', '56'],
    ['10', '56'],
    ['58', '59'],
    ['60', '61'],
    ['58', '61'],
    ['6', '57'],
    ['59', '48'],
    ['5', '60'],
    ['59', '60'],
    ['27', '57'],
    ['57', '28'],
    ['30', '58'],
    ['29', '61'],
    ['5', '61'],
    ['48', '58'],
    ['58', '60'],
    ['59', '61'],
    ['5', '59'],
    ['48', '60'],
    ['29', '58'],
    ['30', '61'],
    ['62', '63'],
    ['73', '70'],
    ['72', '71'],
    ['73', '71'],
    ['72', '70'],
    ['64', '38'],
    ['66', '67'],
    ['68', '69'],
    ['66', '69'],
    ['65', '49'],
    ['62', '68'],
    ['67', '68'],
    ['64', '49'],
    ['63', '69'],
    ['62', '69'],
    ['65', '38'],
    ['66', '68'],
    ['67', '69'],
    ['62', '67'],
    ['49', '38'],
    ['63', '66'],
    ['65', '64'],
    ['39', '62'],
    ['39', '63'],
    ['24', '62'],
    ['5', '62'],
    ['29', '63'],
    ['62', '25'],
    ['62', '29'],
    ['5', '63'],
    ['63', '32'],
    ['3', '63'],
    ['62', '3'],
    ['62', '33'],
    ['67', '39'],
    ['3', '66'],
    ['39', '66'],
    ['39', '68'],
    ['3', '69'],
    ['28', '73'],
    ['28', '70'],
    ['27', '73'],
    ['17', '71'],
    ['37', '72'],
    ['37', '70'],
    ['17', '73'],
    ['77', '75'],
    ['75', '82'],
    ['77', '83'],
    ['74', '81'],
    ['76', '78'],
    ['76', '80'],
    ['74', '79'],
    ['74', '76'],
    ['74', '80'],
    ['76', '79'],
    ['85', '84'],
    ['83', '82'],
    ['85', '82'],
    ['84', '83'],
    ['85', '83'],
    ['84', '82'],
    ['80', '78'],
    ['81', '79'],
    ['80', '79'],
    ['81', '78'],
    ['79', '78'],
    ['81', '80'],
    ['75', '85'],
    ['77', '84'],
    ['77', '82'],
    ['37', '77'],
    ['26', '74'],
    ['34', '77'],
    ['27', '77'],
    ['74', '2'],
    ['3', '74'],
    ['76', '50'],
    ['27', '75'],
    ['32', '76'],
    ['50', '77'],
    ['29', '76'],
    ['76', '34'],
    ['75', '37'],
    ['74', '36'],
    ['30', '74'],
    ['77', '18'],
    ['76', '2'],
    ['74', '50'],
    ['3', '76'],
    ['75', '6'],
    ['6', '77'],
    ['75', '18'],
    ['2', '77'],
    ['57', '75'],
    ['63', '76'],
    ['63', '74'],
    ['18', '84'],
    ['50', '81'],
    ['2', '78'],
    ['2', '81'],
    ['50', '78'],
    ['6', '85'],
    ['6', '82'],
    ['18', '83'],
    ['18', '85'],
    ['53', '54'],
    ['54', '55'],
    ['55', '56'],
    ['53', '56'],
    ['53', '55'],
    ['54', '56'],
  ])),

  dirichlet: {
    '1': bvp.Ux() + bvp.Uz() + bvp.Uy(),
    '4': bvp.Ux() + bvp.Uz() + bvp.Uy(),
    '51': bvp.Ux() + bvp.Uz() + bvp.Uy(),
    '52': bvp.Ux() + bvp.Uz() + bvp.Uy(),
  },

  neumann: {
    '49': bvp.Fz(-5000),
    '59': bvp.Fz(-5000),
    '60': bvp.Fz(-5000),
    '64': bvp.Fz(-5000),
    '67': bvp.Fz(-5000),
    '68': bvp.Fz(-5000),
    '71': bvp.Fz(-5000),
    '72': bvp.Fz(-5000),
    '79': bvp.Fz(-5000),
    '80': bvp.Fz(-5000),
    '83': bvp.Fz(-5000),
    '84': bvp.Fz(-5000),
  },

  expected: {
    tolerance: {
      reaction: 1e-7,
      polynomial: 1e-6,
    },

    reaction: {
      '1': test.Fx(4275.992) + test.Fy(4275.907) + test.Fz(15000),
      '4': test.Fx(-4275.992) + test.Fy(4275.907) + test.Fz(15000),
      '51': test.Fx(-4275.992) + test.Fy(-4275.907) + test.Fz(15000),
      '52': test.Fx(4275.992) + test.Fy(-4275.907) + test.Fz(15000),
    },

    interpolation:
      indexPlusOneToKey(std.map(function(value) test.Constant('Nx', value), [
        -7846.98,
        -0.0858227,
        -7846.956,
        -14095.371,
        4154.644,
        4155.718,
        -14095.371,
        -15003.233,
        1548.565,
        1543.585,
        -14770.659,
        -14805.015,
        -14805.015,
        1311.867,
        -14962.274,
        895.143,
        1301.842,
        -14962.274,
        -13137.605,
        3217.586,
        -4747.775,
        -1244.563,
        -12106.905,
        -1433.553,
        -771.971,
        -1878.868,
        -6867.16,
        -1954.939,
        -2574.812,
        -8310.568,
        -727.126,
        -3037.004,
        -3525.468,
        -3525.468,
        -1219.877,
        -1219.877,
        -443.365,
        -4747.775,
        9229.417,
        -770.904,
        -151.06348,
        -151.06348,
        -14805.015,
        -3753.628,
        1085.305,
        -1334.495,
        -2442.928,
        -2417.909,
        -3840.127,
        -765.006,
        -632.646,
        -632.646,
        -758.467,
        -436.282,
        -1221.361,
        8201.873,
        -14770.659,
        -0.0858227,
        -1244.563,
        -3753.628,
        -151.06348,
        -14770.659,
        -6772.379,
        -7116.386,
        -7116.386,
        -4747.775,
        3217.586,
        -13137.605,
        -14962.274,
        3523.798,
        -7846.956,
        -436.282,
        -4560.937,
        -7846.98,
        3523.798,
        4062.783,
        -7846.98,
        -7846.956,
        -0.0858227,
        -7846.956,
        -14095.371,
        1548.565,
        1543.585,
        -15003.233,
        1085.305,
        1404.208,
        -14805.015,
        4062.783,
        10519.553,
        -443.365,
        -727.126,
        -727.126,
        -770.904,
        8124.475,
        8124.475,
        -5212.373,
        -5212.373,
        -1878.868,
        -771.971,
        -792.749,
        -1221.361,
        -1221.361,
        -4560.937,
        4163.615,
        -0.0858227,
        -14095.371,
        4155.718,
        4154.644,
        -15003.233,
        -14770.659,
        1106.163,
        -14962.274,
        895.143,
        1301.842,
        1400.002,
        -4304.366,
        -12106.905,
        -6867.16,
        -8310.568,
        5981.538,
        -2574.812,
        -3753.628,
        -770.904,
        -3525.468,
        -727.126,
        -443.365,
        10519.553,
        4578.186,
        12144.522,
        1404.208,
        -15003.233,
        -0.0858227,
        -7846.98,
        -5290.423,
        -13137.605,
        -1884.447,
        -1884.447,
        -1334.495,
        -2442.928,
        -2417.909,
        -3840.127,
        -792.749,
        -792.749,
        -765.006,
        -758.467,
        -436.282,
        -1221.361,
        -12106.905,
        -771.971,
        -151.06348,
        1106.163,
        -770.904,
        -443.365,
        1311.867,
        -1219.877,
        -1219.877,
        -6739.727,
        -3037.004,
        -6739.727,
        7510.323,
        7510.323,
        -758.467,
        -758.467,
        -0.0858227,
        1400.002,
        -1334.495,
        -765.006,
        -765.006,
        -632.646,
        -632.646,
        -436.282,
        -0.0858227,
        -4561.114,
        -4561.114,
        -4561.114,
        -4561.114,
        -0.0858227,
        -4560.937,
        -4560.937,
        -4545.782,
        -4545.782,
        -4545.605,
        -4545.605,
        -4545.782,
        -4545.782,
        -4545.605,
        -4545.605,
        7222.761,
        7222.761,
        -1057.925,
        0,
        -6817.052,
        -6817.052,
        -1396.304,
        0,
        0,
        8176.207,
        8176.207,
        3518.728,
        3518.728,
        5618.292,
        5618.292,
        -5426.819,
        -5426.819,
        4346.917,
        4346.917,
        -12106.905,
        -591.397,
        -1476.169,
        5888.002,
        5888.002,
        5492.022,
        7222.761,
        7222.761,
        -1057.925,
        5492.022,
        -6817.052,
        -1396.304,
        -1422.727,
        8176.207,
        3518.728,
        -842.585,
        5618.292,
        5618.292,
        -5426.819,
        7044.708,
        4346.917,
        7044.708,
        4163.615,
        -3840.127,
        -13137.605,
        -5290.423,
        12144.522,
        -3753.628,
        -2574.812,
        -2574.812,
        -1954.939,
        -1433.553,
        -3840.127,
        -792.749,
        -6817.052,
        8176.207,
        3518.728,
        -5426.819,
        4346.917,
        8201.873,
        4130.187,
        4130.187,
        -5644.96,
        -5644.96,
        4710.306,
        4710.306,
        -4747.775,
        8201.873,
        -7116.386,
        3523.798,
        3523.798,
        -5212.373,
        -5212.373,
        5981.538,
        -6739.727,
        -6739.727,
        7510.323,
        7510.323,
        -591.397,
        -1476.169,
        5888.002,
        5888.002,
        5492.022,
        5492.022,
        -1422.727,
        -842.585,
        7044.708,
        7044.708,
        4130.187,
        -5644.96,
        4710.306,
        -6772.379,
        -5694.159,
        -771.971,
        -1244.563,
        -8310.568,
        -6867.16,
        -8310.568,
        9229.417,
        -5694.159,
        -3525.468,
        -1954.939,
        -1878.868,
        -1244.563,
        -1878.868,
        -1954.939,
        4578.186,
        -2442.928,
        -2442.928,
        -2417.909,
        -4304.366,
        -1884.447,
        -1884.447,
        -1334.495,
        0,
        -6867.16,
        -2417.909,
        -7116.386,
        4062.783,
        4062.783,
        8124.475,
        8124.475,
        8201.873,
        4130.187,
        -5644.96,
        4710.306,
        -0.033728,
        -0.033728,
        -0.033728,
        -0.033728,
        -0.033729,
        -0.033727,
      ])),

  },
}
