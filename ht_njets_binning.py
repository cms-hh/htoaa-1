from collections import OrderedDict as OD

htbin = OD([
    ('HT-0To70', [0, 70]),
    ('HT-70To100', [70,100]),
    ('HT-100To200', [100,200]),
    ('HT-200To400', [200,400]),
    ('HT-400To600', [400,600]),
    ('HT-600To800', [600,800]),
    ('HT-800To1200', [800, 1200]),
    ('HT-1200To2500', [1200,2500]),
    ('HT-2500ToInf', [2500]),
])
njetbin = OD([
    ('0Jet', 0),
    ('1Jet', 1),
    ('2Jet', 2),
    ('3Jet', 3),
    ('4Jet', 4)
])
