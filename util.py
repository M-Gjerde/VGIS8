import numpy as np

arr = [(-0.30920000000000003, 0.5078), (-15.0655, 0.6555000000000001), (-8.5904, 0.4539), (-0.25, 0),
         (-5.9583, 0.5383), (-7.5986, 0.40030000000000004), (-14.6286, 0.38120000000000004),
         (-2.5068, 0.41100000000000003), (-13.0997, 0.7562), (-14.834100000000001, 0.5666), (-0.2194, 8.6981),
         (-19.6679, 1.5250000000000001), (-0.13390000000000002, 0), (-0.37270000000000003, 0.3546),
         (-0.2884, 0.44520000000000004), (-3.0599000000000003, 0.6312), (-9.0771, 0.5111),
         (-0.6938000000000001, 0.36100000000000004), (-19.5383, 0.5511), (-0.3315, 0.5249),
         (-0.1557, 0.5579000000000001), (-29.8609, 7.944400000000001), (-4.0993, 0),
         (-0.6541, 0.36920000000000003), (-2.4477, 0.562), (-0.46340000000000003, 0.4247), (-0.393, 0.4307),
         (-1.3047, 0.3443), (-0.47150000000000003, 0), (-0.392, 1.9998),
         (-6.594200000000001, 0), (-4.661300000000001, 152.9713), (-11.3763, 2.5164),
         (-0.509, 87.5239), (-0.3972, 0), (-17.346700000000002, 0.3705), (-12.6111, 0.4546),
         (-8.3622, 0), (-10.455400000000001, 0.7605000000000001), (-2.0926, 0.5157),
         (-0.47350000000000003, 0.9341), (-1.8532000000000002, 0.453), (-11.740400000000001, 0.369),
         (-26.8673, 0.6022000000000001), (-0.21680000000000002, 0.5725), (-0.45730000000000004, 0.342),
         (-12.4122, 0.3593), (-0.3766, 0.2504), (-24.6383, 0.5058), (-14.5701, 1.8908), (-0.376, 5.7869),
         (-12.611600000000001, 1.5661), (-0.22, 0.39380000000000004), (-0.2482, 0.33440000000000003),
         (-6.2967, 0.44620000000000004), (-6.9081, 0.2561), (-10.395900000000001, 6.0624), (-12.394, 0.2792),
         (-4.1367, 0.375), (-0.5882000000000001, 2.4867), (-0.0279, 0.11570000000000001)]

arr = [(-41.5422, -41.0893), (-46.075, -29.2437), (-35.8853, -34.2588), (-46.874900000000004, -32.8872), (-42.294200000000004, -34.214600000000004), (-40.387100000000004, -36.134), (-47.4846, -23.740000000000002), (-41.543, -41.0934), (-46.0713, -29.2364), (-35.8985, -34.272200000000005), (-46.8746, -32.8933), (-42.3012, -34.2381), (-40.3952, -36.1306), (-47.4851, -23.7498), (-41.5403, -41.0908), (-46.0715, -29.2469), (-35.8853, -34.2671), (-46.886500000000005, -32.9032), (-42.299, -34.22), (-40.3956, -36.1186), (-47.481700000000004, -23.7236), (-41.5403, -41.077200000000005), (-38.937400000000004, -29.161800000000003), (-46.0739, -29.1342), (-35.873200000000004, -34.2541), (-46.897400000000005, -32.890100000000004)]

np.array(arr)

print(np.mean(arr, axis=0))