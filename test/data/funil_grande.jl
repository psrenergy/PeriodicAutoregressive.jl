# Funil grande from 1931 to 2019
const funil_grande = [302.0, 537.0, 298.0, 167.0, 80.0, 43.4, 85.0, 76.0, 75.0, 35.0, 76.0, 158.0,
    159.0, 261.0, 282.0, 220.0, 163.0, 105.0, 84.0, 71.0, 73.0, 96.0, 110.0, 327.0,
    407.0, 198.0, 178.0, 119.0, 86.0, 72.0, 68.0, 62.0, 59.0, 90.0, 91.0, 169.0,
    278.0, 129.0, 141.0, 86.0, 67.0, 58.0, 53.0, 46.0, 51.0, 57.0, 81.0, 252.0,
    383.0, 539.0, 263.0, 227.0, 158.0, 115.0, 106.0, 81.0, 69.0, 100.0, 100.0, 117.0,
    81.0, 102.0, 259.0, 163.0, 92.0, 67.0, 58.0, 59.0, 62.0, 58.0, 129.0, 186.0,
    540.0, 333.0, 227.0, 155.0, 137.0, 98.0, 78.0, 63.0, 57.0, 163.0, 253.0, 561.0,
    367.0, 397.0, 329.0, 227.0, 182.0, 129.0, 100.0, 91.0, 89.0, 144.0, 199.0, 404.0,
    466.0, 393.0, 231.0, 207.0, 138.0, 101.0, 85.0, 67.0, 60.0, 81.0, 82.0, 256.0,
    321.0, 432.0, 340.0, 214.0, 139.0, 102.0, 81.0, 65.0, 66.0, 79.0, 217.0, 269.0,
    444.0, 270.0, 267.0, 208.0, 135.0, 116.0, 110.0, 85.0, 107.0, 124.0, 158.0, 344.0,
    327.0, 263.0, 372.0, 207.0, 160.0, 122.0, 101.0, 86.0, 80.0, 118.0, 182.0, 341.0,
    734.0, 556.0, 511.0, 320.0, 204.0, 111.0, 96.0, 79.0, 74.0, 96.0, 91.0, 216.0,
    194.0, 342.0, 327.0, 180.0, 119.0, 84.0, 72.0, 58.0, 50.0, 55.0, 71.0, 128.0,
    184.0, 286.0, 191.0, 177.0, 102.0, 83.0, 68.0, 53.0, 51.0, 55.0, 93.0, 243.0,
    460.0, 180.0, 209.0, 164.0, 110.0, 87.0, 73.0, 63.0, 56.0, 82.0, 148.0, 143.0,
    270.0, 261.0, 686.0, 286.0, 186.0, 142.0, 115.0, 100.0, 119.0, 101.0, 108.0, 268.0,
    270.0, 280.0, 298.0, 181.0, 125.0, 99.0, 79.0, 65.0, 58.0, 57.0, 111.0, 281.0,
    298.0, 399.0, 237.0, 167.0, 113.0, 120.0, 83.0, 64.0, 56.0, 75.0, 76.0, 231.0,
    269.0, 392.0, 240.0, 189.0, 141.0, 106.0, 86.0, 70.0, 65.0, 83.0, 261.0, 327.0,
    294.0, 447.0, 396.0, 273.0, 168.0, 135.0, 110.0, 94.0, 79.0, 88.0, 74.0, 160.0,
    274.0, 369.0, 415.0, 259.0, 158.0, 136.0, 106.0, 85.0, 71.0, 76.0, 138.0, 212.0,
    152.0, 149.0, 145.0, 185.0, 113.0, 89.0, 69.0, 46.0, 61.0, 61.0, 98.0, 170.0,
    109.0, 173.0, 101.0, 125.0, 84.0, 45.0, 54.0, 46.0, 41.0, 54.0, 90.0, 127.0,
    172.0, 129.0, 137.0, 110.0, 71.0, 63.0, 40.0, 33.1, 27.7, 33.0, 47.0, 158.0,
    188.0, 135.0, 148.0, 107.0, 88.0, 78.0, 60.0, 61.0, 53.0, 33.0, 60.0, 244.0,
    218.0, 193.0, 213.0, 249.0, 137.0, 117.0, 93.0, 74.0, 72.0, 54.0, 97.0, 220.0,
    117.0, 115.0, 153.0, 115.0, 121.0, 82.0, 66.0, 50.0, 56.0, 44.0, 40.0, 90.0,
    257.0, 201.0, 142.0, 146.0, 126.0, 112.0, 98.0, 79.0, 70.0, 84.0, 100.0, 131.0,
    285.0, 321.0, 404.0, 177.0, 158.0, 120.0, 135.0, 106.0, 76.0, 81.0, 95.0, 225.0,
    453.0, 497.0, 480.0, 265.0, 234.0, 160.0, 132.0, 114.0, 106.0, 85.0, 125.0, 144.0,
    256.0, 361.0, 311.0, 181.0, 161.0, 153.0, 117.0, 103.0, 92.0, 146.0, 159.0, 231.0,
    329.0, 327.0, 188.0, 124.0, 97.0, 78.0, 71.0, 59.0, 52.0, 73.0, 100.0, 74.1,
    316.0, 392.0, 145.0, 172.0, 135.0, 103.0, 96.0, 80.0, 68.0, 105.0, 145.0, 246.0,
    509.0, 671.0, 499.0, 336.0, 286.0, 217.0, 192.0, 158.0, 179.0, 258.0, 313.0, 344.0,
    607.0, 310.0, 362.0, 228.0, 148.0, 119.0, 103.0, 88.0, 79.0, 185.0, 270.0, 367.0,
    550.0, 416.0, 285.0, 186.0, 146.0, 130.0, 106.0, 97.0, 78.0, 92.0, 231.0, 244.0,
    298.0, 189.0, 191.0, 124.0, 98.0, 81.0, 74.0, 78.0, 90.0, 111.0, 105.0, 311.0,
    270.0, 255.0, 208.0, 148.0, 98.0, 108.0, 89.0, 81.0, 74.0, 135.0, 272.0, 259.0,
    300.0, 207.0, 172.0, 144.0, 95.0, 79.0, 80.0, 74.0, 85.0, 121.0, 191.0, 131.0,
    150.0, 96.0, 102.0, 74.0, 65.0, 82.0, 67.0, 63.0, 67.0, 96.0, 168.0, 373.0,
    211.0, 290.0, 310.0, 174.0, 129.0, 102.0, 118.0, 92.0, 86.0, 126.0, 300.0, 280.0,
    341.0, 332.0, 155.0, 195.0, 148.0, 110.0, 106.0, 100.0, 97.0, 115.0, 170.0, 257.0,
    320.0, 228.0, 283.0, 238.0, 147.0, 111.0, 96.0, 82.0, 66.0, 102.0, 102.0, 208.0,
    294.0, 265.0, 154.0, 115.0, 98.0, 81.0, 84.0, 70.0, 60.0, 98.0, 214.0, 257.0,
    174.0, 184.0, 205.0, 138.0, 122.0, 110.0, 114.0, 118.0, 169.0, 171.0, 226.0, 330.0,
    283.0, 240.0, 223.0, 192.0, 121.0, 109.0, 94.0, 80.0, 103.0, 85.0, 160.0, 251.0,
    495.0, 226.0, 200.0, 156.0, 127.0, 120.0, 104.0, 86.0, 85.0, 110.0, 187.0, 221.0,
    323.0, 684.0, 329.0, 184.0, 160.0, 134.0, 127.0, 109.0, 137.0, 120.0, 190.0, 338.0,
    509.0, 321.0, 191.0, 238.0, 144.0, 126.0, 111.0, 96.0, 89.0, 95.0, 152.0, 351.0,
    420.0, 229.0, 223.0, 144.0, 119.0, 113.0, 88.0, 92.0, 81.0, 120.0, 276.0, 352.0,
    402.0, 259.0, 424.0, 268.0, 167.0, 142.0, 123.0, 113.0, 91.0, 160.0, 206.0, 425.0,
    625.0, 462.0, 508.0, 377.0, 222.0, 217.0, 182.0, 151.0, 181.0, 283.0, 341.0, 551.0,
    318.0, 194.0, 180.0, 145.0, 140.0, 102.0, 81.0, 71.0, 96.0, 94.0, 129.0, 303.0,
    515.0, 422.0, 434.0, 219.0, 155.0, 128.0, 108.0, 90.0, 95.0, 103.0, 163.0, 221.0,
    354.0, 305.0, 225.0, 139.0, 122.0, 92.0, 90.0, 100.0, 79.0, 62.0, 74.0, 353.0,
    349.0, 310.0, 231.0, 175.0, 143.0, 120.0, 99.0, 75.0, 97.0, 88.0, 112.0, 315.0,
    301.0, 416.0, 249.0, 168.0, 133.0, 114.0, 87.0, 75.0, 63.0, 116.0, 160.0, 206.0,
    341.0, 305.0, 308.0, 174.0, 123.0, 116.0, 97.0, 93.0, 94.0, 75.0, 89.0, 205.0,
    251.0, 128.0, 199.0, 146.0, 119.0, 91.0, 89.0, 79.0, 97.0, 94.0, 111.0, 134.0,
    485.0, 367.0, 296.0, 260.0, 153.0, 121.0, 100.0, 81.0, 72.0, 118.0, 123.0, 172.0,
    803.0, 467.0, 231.0, 171.0, 145.0, 100.0, 82.0, 73.0, 95.0, 110.0, 231.0, 182.0,
    234.0, 302.0, 289.0, 217.0, 133.0, 126.0, 92.0, 78.0, 77.0, 123.0, 103.0, 159.0,
    451.0, 208.0, 282.0, 183.0, 191.0, 130.0, 103.0, 83.0, 69.0, 79.0, 95.0, 216.0,
    172.0, 400.0, 215.0, 161.0, 129.0, 103.0, 87.0, 65.0, 63.0, 91.0, 130.0, 215.0,
    349.0, 249.0, 241.0, 162.0, 122.0, 106.0, 95.0, 74.0, 99.0, 96.0, 312.0, 360.0,
    839.0, 318.0, 305.0, 188.0, 137.0, 130.0, 104.0, 84.0, 77.0, 85.0, 123.0, 210.0,
    230.0, 235.0, 179.0, 133.0, 105.0, 96.0, 75.0, 76.0, 59.0, 86.0, 133.0, 212.0,
    264.0, 211.0, 279.0, 151.0, 105.0, 91.0, 68.0, 54.0, 45.0, 49.0, 78.0, 162.0,
    329.0, 305.0, 224.0, 152.0, 101.0, 89.0, 77.0, 66.0, 90.0, 63.0, 129.0, 171.0,
    172.0, 125.0, 123.0, 116.0, 75.0, 61.0, 52.0, 49.0, 54.0, 67.0, 129.0, 241.0,
    246.0, 348.0, 205.0, 116.0, 87.0, 96.0, 84.0, 67.0, 62.0, 43.0, 87.0, 134.0,
    329.0, 223.0, 193.0, 114.0, 88.0, 67.0, 58.0, 54.0, 51.0, 45.0, 75.0, 198.0,
    263.0, 296.0, 254.0, 208.0, 123.0, 110.0, 94.0, 70.0, 73.0, 62.0, 84.0, 307.0,
    428.0, 273.0, 278.0, 161.0, 145.0, 115.0, 93.0, 74.0, 75.0, 64.0, 149.0, 341.0,
    190.0, 222.0, 220.0, 124.0, 102.0, 83.0, 71.0, 58.0, 64.0, 96.0, 130.0, 215.0,
    498.0, 307.0, 162.0, 124.0, 98.0, 86.0, 73.0, 56.0, 43.0, 54.0, 92.0, 160.0,
    197.0, 308.0, 339.0, 232.0, 128.0, 103.0, 82.0, 71.0, 81.0, 90.0, 159.0, 391.0,
    454.0, 355.0, 299.0, 304.0, 152.0, 127.0, 104.0, 79.0, 97.0, 160.0, 151.0, 351.0,
    345.0, 191.0, 298.0, 172.0, 124.0, 102.0, 88.0, 67.0, 61.0, 98.0, 211.0, 326.0,
    430.0, 165.0, 321.0, 197.0, 120.0, 112.0, 87.0, 68.0, 50.0, 79.0, 100.0, 417.0,
    685.0, 274.0, 222.0, 188.0, 149.0, 138.0, 108.0, 85.0, 69.0, 68.0, 102.0, 126.0,
    270.0, 250.0, 207.0, 178.0, 113.0, 111.0, 83.0, 62.0, 63.0, 85.0, 95.0, 213.0,
    123.0, 72.0, 79.0, 86.0, 56.0, 50.0, 44.0, 41.0, 31.0, 24.0, 88.0, 111.0,
    75.4, 136.0, 187.0, 118.0, 83.0, 64.0, 50.0, 40.0, 75.0, 38.0, 118.0, 201.0,
    288.0, 210.0, 187.0, 96.0, 78.0, 82.0, 58.0, 48.0, 46.0, 54.0, 141.0, 174.0,
    155.0, 115.0, 106.0, 74.0, 74.0, 63.0, 47.0, 37.0, 29.0, 46.0, 60.0, 144.0,
    170.0, 141.0, 210.0, 88.0, 59.0, 55.0, 44.0, 57.0, 43.0, 73.0, 135.0, 197.0,
    134.0, 145.0, 215.0, 127.0, 87.0, 69.0, 54.0, 47.0, 40.0, 45.0, 100.0, 158.0]
