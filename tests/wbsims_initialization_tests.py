import unittest

from wbsims_init import coalsims_fx



gt_array, mutations = coalsims_fx([10],1,.0001,[0],[5],13000,7.6E-8,0,1800,23,240)
gt_array, mutations = coalsims_fx([10,10],2,.0001,[1000],[5,5],13000,7.6E-8,0,1800,23,240)
#gt_array, gt_array2, mutations = coalsims_fx([10],1,.0001,[0],[5],13000,7.6E-8,2.9E-9,1800,23,240)
#gt_array, gt_array2, mutations = coalsims_fx([10,10],2,.0001,[1000],[5,5],13000,7.6E-8,2.9E-9,1800,23,240)
