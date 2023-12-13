CDF  w   
      time       lon       lat          	   CDI       @Climate Data Interface version 2.0.5 (https://mpimet.mpg.de/cdi)   Conventions       CF-1.5     source        #PISM development no-version-control    command      / /home/m/m300792/ModelCodes/mPISM_build/install_03-10-22/bin/pismr -i pism_-065000/pism_-065000.nc -stress_balance ssa+sia -y 100 -ys -65000 -o pism_-064900/pism_-064900.nc -surface_given_file forcing/smb_-064900.nc -surface given -ocean_th_file forcing/oce_-064900.nc -ocean th -topg_file pism_-065000/pism_-065000_new_topg.nc -bed_def prescribed -uplift_file pism_-065000/pism_-065000_dbdt.nc -extra_times 10 -extra_vars thk,usurf,bheatflx,velsurf_mag,velbase_mag,tillwat,bmelt,velbase,flux_divergence,velsurf,velbar,wvelsurf,wvelbase,climatic_mass_balance_cumulative,potential_climatic_mass_balance_cumulative,tempbase,tempsurf,mask,temppabase,tempicethk_basal,topg,nonneg_flux_cumulative,grounded_basal_flux_cumulative,floating_basal_flux_cumulative,discharge_flux_cumulative,Href,taud_mag,lat,lon,calving_mask,sea_level -extra_file pism_-064900/pism_-064900_extra -extra_split -ts_times yearly -ts_file pism_-064900/pism_-064900_ts.nc -topg_to_phi 15.0,30.0,-300.0,400.0 -ssafd_ksp_rtol 1e-7 -enthalpy_temperate_conductivity_ratio 0.001 -pseudo_plastic -pseudo_plastic_q .25 -pseudo_plastic_uthreshold 70 -till_effective_fraction_overburden 0.04 -tauc_slippery_grounding_lines -sia_e 8 -ssa_e 1 -ssa_upwind_factor 0.7 -calving eigen_calving,thickness_calving -thickness_calving_threshold 100 -part_redist -part_grid -cfbc -kill_icebergs -eigen_calving_K 1e16 -verbose 2 -sliding_scale_factor_reduces_tauc 1 -th_heat_coefficient 4e-5 -th_salinity_coefficient 4e-8 -o_format netcdf3 -sediment_file till_mask.nc -no_divide_tillphi -periodicity none -ocean_th_period 1 -o_order yxz
    proj4         j+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs     history_of_appended_files        �Tue Oct 17 11:11:31 2023: Appended file outFinal2 had following "history" attribute:
Tue Oct 17 11:11:30 2023: ncap2 -O -s lon_bnds=float(lon_bnds) outFinal1 outFinal2
Tue Oct 17 11:11:29 2023: ncap2 -O -s lat_bnds=float(lat_bnds) outFinal outFinal1
Tue Oct 17 11:11:29 2023: cdo -add thkNew thkOld outFinal
Tue Oct 17 11:11:27 2023: cdo -mul -selvar,thk pism_-064900/pism_-064900.nc MaskLaurentide thkNew
     NCO       _netCDF Operators version 5.0.6 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)    history      1�Sun Nov 05 20:53:38 2023: cdo -mergetime File000001.nc File000002.nc File000003.nc File000004.nc File000005.nc File000006.nc File000007.nc File000008.nc File000009.nc File000010.nc File000011.nc File000012.nc File000013.nc File000014.nc File000015.nc File000016.nc File000017.nc File000018.nc File000019.nc File000020.nc File000021.nc File000022.nc File000023.nc File000024.nc File000025.nc File000026.nc File000027.nc File000028.nc File000029.nc File000030.nc File000031.nc File000032.nc File000033.nc File000034.nc File000035.nc File000036.nc File000037.nc File000038.nc File000039.nc File000040.nc File000041.nc File000042.nc File000043.nc File000044.nc File000045.nc File000046.nc File000047.nc File000048.nc File000049.nc File000050.nc File000051.nc File000052.nc File000053.nc File000054.nc File000055.nc File000056.nc File000057.nc File000058.nc File000059.nc File000060.nc File000061.nc File000062.nc File000063.nc File000064.nc File000065.nc File000066.nc File000067.nc File000068.nc File000069.nc File000070.nc File000071.nc File000072.nc File000073.nc File000074.nc File000075.nc File000076.nc File000077.nc File000078.nc File000079.nc File000080.nc File000081.nc File000082.nc File000083.nc File000084.nc File000085.nc File000086.nc File000087.nc File000088.nc File000089.nc File000090.nc File000091.nc File000092.nc File000093.nc File000094.nc File000095.nc File000096.nc File000097.nc File000098.nc File000099.nc File000100.nc File000101.nc File000102.nc File000103.nc File000104.nc File000105.nc File000106.nc File000107.nc File000108.nc File000109.nc File000110.nc File000111.nc File000112.nc File000113.nc File000114.nc File000115.nc File000116.nc File000117.nc File000118.nc File000119.nc File000120.nc File000121.nc File000122.nc File000123.nc File000124.nc File000125.nc File000126.nc File000127.nc File000128.nc File000129.nc File000130.nc File000131.nc File000132.nc File000133.nc File000134.nc File000135.nc File000136.nc File000137.nc File000138.nc File000139.nc File000140.nc File000141.nc File000142.nc File000143.nc File000144.nc File000145.nc File000146.nc File000147.nc File000148.nc File000149.nc File000150.nc File000151.nc File000152.nc File000153.nc File000154.nc File000155.nc File000156.nc File000157.nc File000158.nc File000159.nc File000160.nc File000161.nc File000162.nc File000163.nc File000164.nc File000165.nc File000166.nc File000167.nc File000168.nc File000169.nc File000170.nc File000171.nc File000172.nc File000173.nc File000174.nc File000175.nc File000176.nc File000177.nc File000178.nc File000179.nc File000180.nc File000181.nc File000182.nc File000183.nc File000184.nc File000185.nc File000186.nc File000187.nc File000188.nc File000189.nc File000190.nc File000191.nc File000192.nc File000193.nc File000194.nc File000195.nc File000196.nc File000197.nc File000198.nc File000199.nc File000200.nc File000201.nc File000202.nc File000203.nc File000204.nc File000205.nc File000206.nc File000207.nc File000208.nc File000209.nc File000210.nc File000211.nc File000212.nc File000213.nc File000214.nc File000215.nc File000216.nc File000217.nc File000218.nc File000219.nc File000220.nc File000221.nc File000222.nc File000223.nc File000224.nc File000225.nc File000226.nc File000227.nc File000228.nc File000229.nc File000230.nc File000231.nc File000232.nc File000233.nc File000234.nc File000235.nc File000236.nc File000237.nc File000238.nc File000239.nc File000240.nc File000241.nc File000242.nc File000243.nc File000244.nc File000245.nc File000246.nc File000247.nc File000248.nc File000249.nc File000250.nc File000251.nc File000252.nc File000253.nc File000254.nc File000255.nc File000256.nc File000257.nc File000258.nc File000259.nc File000260.nc File000261.nc File000262.nc File000263.nc File000264.nc File000265.nc File000266.nc File000267.nc File000268.nc File000269.nc File000270.nc File000271.nc File000272.nc File000273.nc File000274.nc File000275.nc File000276.nc File000277.nc File000278.nc File000279.nc File000280.nc File000281.nc File000282.nc File000283.nc File000284.nc File000285.nc File000286.nc File000287.nc File000288.nc File000289.nc File000290.nc File000291.nc File000292.nc File000293.nc File000294.nc File000295.nc File000296.nc File000297.nc File000298.nc File000299.nc File000300.nc File000301.nc File000302.nc File000303.nc File000304.nc File000305.nc File000306.nc File000307.nc File000308.nc File000309.nc File000310.nc File000311.nc File000312.nc File000313.nc File000314.nc File000315.nc File000316.nc File000317.nc File000318.nc File000319.nc File000320.nc File000321.nc File000322.nc File000323.nc File000324.nc File000325.nc File000326.nc File000327.nc File000328.nc File000329.nc File000330.nc File000331.nc File000332.nc File000333.nc File000334.nc File000335.nc File000336.nc File000337.nc File000338.nc File000339.nc File000340.nc File000341.nc File000342.nc File000343.nc File000344.nc File000345.nc File000346.nc File000347.nc File000348.nc File000349.nc File000350.nc File000351.nc File000352.nc File000353.nc File000354.nc File000355.nc File000356.nc File000357.nc File000358.nc File000359.nc File000360.nc File000361.nc File000362.nc File000363.nc File000364.nc File000365.nc File000366.nc File000367.nc File000368.nc File000369.nc File000370.nc File000371.nc File000372.nc File000373.nc File000374.nc File000375.nc File000376.nc File000377.nc File000378.nc File000379.nc File000380.nc File000381.nc File000382.nc File000383.nc File000384.nc File000385.nc File000386.nc File000387.nc File000388.nc File000389.nc File000390.nc File000391.nc File000392.nc File000393.nc File000394.nc File000395.nc File000396.nc File000397.nc File000398.nc File000399.nc File000400.nc File000401.nc File000402.nc File000403.nc File000404.nc File000405.nc File000406.nc File000407.nc File000408.nc File000409.nc File000410.nc File000411.nc File000412.nc File000413.nc File000414.nc File000415.nc File000416.nc File000417.nc File000418.nc File000419.nc File000420.nc File000421.nc File000422.nc File000423.nc File000424.nc File000425.nc File000426.nc File000427.nc File000428.nc File000429.nc File000430.nc File000431.nc File000432.nc File000433.nc File000434.nc File000435.nc File000436.nc File000437.nc File000438.nc File000439.nc File000440.nc File000441.nc File000442.nc File000443.nc File000444.nc File000445.nc File000446.nc File000447.nc File000448.nc File000449.nc File000450.nc File000451.nc File000452.nc File000453.nc File000454.nc File000455.nc File000456.nc File000457.nc File000458.nc File000459.nc File000460.nc File000461.nc File000462.nc File000463.nc File000464.nc File000465.nc File000466.nc File000467.nc File000468.nc File000469.nc File000470.nc File000471.nc File000472.nc File000473.nc File000474.nc File000475.nc File000476.nc File000477.nc File000478.nc File000479.nc File000480.nc File000481.nc File000482.nc File000483.nc File000484.nc File000485.nc File000486.nc File000487.nc File000488.nc File000489.nc File000490.nc File000491.nc File000492.nc File000493.nc File000494.nc File000495.nc File000496.nc File000497.nc File000498.nc File000499.nc File000500.nc File000501.nc File000502.nc File000503.nc File000504.nc File000505.nc File000506.nc File000507.nc File000508.nc File000509.nc File000510.nc File000511.nc File000512.nc File000513.nc File000514.nc File000515.nc File000516.nc File000517.nc File000518.nc File000519.nc File000520.nc File000521.nc File000522.nc File000523.nc File000524.nc File000525.nc File000526.nc File000527.nc File000528.nc File000529.nc File000530.nc File000531.nc File000532.nc File000533.nc File000534.nc File000535.nc File000536.nc File000537.nc File000538.nc File000539.nc File000540.nc File000541.nc File000542.nc File000543.nc File000544.nc File000545.nc File000546.nc File000547.nc File000548.nc File000549.nc File000550.nc File000551.nc File000552.nc File000553.nc File000554.nc File000555.nc File000556.nc File000557.nc File000558.nc File000559.nc File000560.nc File000561.nc File000562.nc File000563.nc File000564.nc File000565.nc File000566.nc File000567.nc File000568.nc File000569.nc File000570.nc File000571.nc File000572.nc File000573.nc File000574.nc File000575.nc File000576.nc File000577.nc File000578.nc File000579.nc File000580.nc File000581.nc File000582.nc File000583.nc File000584.nc File000585.nc File000586.nc File000587.nc File000588.nc File000589.nc File000590.nc File000591.nc File000592.nc File000593.nc File000594.nc File000595.nc File000596.nc File000597.nc File000598.nc File000599.nc File000600.nc File000601.nc File000602.nc File000603.nc File000604.nc File000605.nc File000606.nc File000607.nc File000608.nc File000609.nc File000610.nc File000611.nc File000612.nc File000613.nc File000614.nc File000615.nc File000616.nc File000617.nc File000618.nc File000619.nc File000620.nc File000621.nc File000622.nc File000623.nc File000624.nc File000625.nc File000626.nc File000627.nc File000628.nc File000629.nc File000630.nc File000631.nc File000632.nc File000633.nc File000634.nc File000635.nc File000636.nc File000637.nc File000638.nc File000639.nc File000640.nc File000641.nc File000642.nc File000643.nc File000644.nc File000645.nc File000646.nc File000647.nc File000648.nc File000649.nc File000650.nc File000651.nc File000652.nc File000653.nc File000654.nc File000655.nc File000656.nc File000657.nc File000658.nc File000659.nc File000660.nc File000661.nc File000662.nc File000663.nc File000664.nc File000665.nc File000666.nc File000667.nc File000668.nc File000669.nc File000670.nc File000671.nc File000672.nc File000673.nc File000674.nc File000675.nc File000676.nc File000677.nc File000678.nc File000679.nc File000680.nc File000681.nc File000682.nc File000683.nc File000684.nc File000685.nc File000686.nc File000687.nc File000688.nc File000689.nc File000690.nc File000691.nc File000692.nc File000693.nc File000694.nc File000695.nc File000696.nc File000697.nc File000698.nc File000699.nc File000700.nc File000701.nc File000702.nc File000703.nc File000704.nc File000705.nc File000706.nc File000707.nc File000708.nc File000709.nc File000710.nc File000711.nc File000712.nc File000713.nc File000714.nc File000715.nc File000716.nc File000717.nc File000718.nc File000719.nc File000720.nc File000721.nc File000722.nc File000723.nc File000724.nc File000725.nc File000726.nc File000727.nc File000728.nc File000729.nc File000730.nc File000731.nc File000732.nc File000733.nc File000734.nc File000735.nc File000736.nc File000737.nc File000738.nc File000739.nc File000740.nc File000741.nc File000742.nc File000743.nc File000744.nc File000745.nc File000746.nc File000747.nc File000748.nc File000749.nc File000750.nc File000751.nc File000752.nc File000753.nc File000754.nc File000755.nc File000756.nc File000757.nc File000758.nc File000759.nc File000760.nc File000761.nc File000762.nc File000763.nc File000764.nc File000765.nc File000766.nc File000767.nc File000768.nc File000769.nc File000770.nc File000771.nc File000772.nc File000773.nc File000774.nc File000775.nc File000776.nc File000777.nc File000778.nc File000779.nc File000780.nc File000781.nc File000782.nc File000783.nc File000784.nc File000785.nc File000786.nc File000787.nc File000788.nc File000789.nc File000790.nc File000791.nc File000792.nc File000793.nc File000794.nc File000795.nc File000796.nc File000797.nc File000798.nc File000799.nc File000800.nc File000801.nc File000802.nc File000803.nc File000804.nc File000805.nc File000806.nc File000807.nc File000808.nc File000809.nc File000810.nc File000811.nc File000812.nc File000813.nc File000814.nc File000815.nc File000816.nc File000817.nc File000818.nc File000819.nc File000820.nc File000821.nc File000822.nc File000823.nc File000824.nc File000825.nc File000826.nc File000827.nc File000828.nc File000829.nc File000830.nc File000831.nc File000832.nc File000833.nc File000834.nc File000835.nc File000836.nc File000837.nc File000838.nc File000839.nc File000840.nc File000841.nc File000842.nc File000843.nc File000844.nc File000845.nc File000846.nc File000847.nc File000848.nc File000849.nc File000850.nc File000851.nc File000852.nc File000853.nc File000854.nc File000855.nc File000856.nc File000857.nc File000858.nc File000859.nc File000860.nc File000861.nc File000862.nc File000863.nc File000864.nc File000865.nc File000866.nc File000867.nc File000868.nc File000869.nc File000870.nc File000871.nc File000872.nc File000873.nc File000874.nc File000875.nc File000876.nc File000877.nc File000878.nc File000879.nc File000880.nc File000881.nc File000882.nc File000883.nc File000884.nc File000885.nc File000886.nc File000887.nc IceVolumeHudson.nc
Sun Nov 05 20:36:37 2023: cdo -fldsum -mulc,10000 -mulc,10000 -selindexbox,279,477,230,362 -selvar,thk /work/ba0989/m300792/HE_Runs/ExperimentsComposite//HE112//pism_-064900/pism_-064900.nc Tmp/File000001.nc      CDO       @Climate Data Operators version 2.0.5 (https://mpimet.mpg.de/cdo)         time                standard_name         time   	long_name         time   units         seconds since 1-1-1    calendar      365_day    axis      T               >�   lon                standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X               >�   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y               >�   thk                       standard_name         land_ice_thickness     	long_name         land ice thickness     units         m      pism_intent       model_state             >�                �}Ȁ@�  C4�(w-o�}���   C4��U���}�=   C4���[��}�A�@  C4����@�}��9`  C4���7��}�·�  C3P;�8[�}�5�  C2{�蒨6�}vC��  C1��!�N��}j�1�  C1` 8��&�}^İ   C1ӛk��}S.   C0�M7, "�}GE�@  C0���j���};�*`  C0�����}/ƨ�  C0u�XX�}$&�  C0m�j�Tw�}G��  C0���1�}�"�  C0��େ$�} ȡ   C0�e�,[��|�	   C0�jط(8�|�I�@  C1�MRmQ�|݊`  C1;��)�>�|�ʙ�  C1Yey���|��  C1sK���|�K��  C1�Pڒ1�|���  C1�-���|�̒   C1� q�	��|�   C1���x��|�M�@  C1ʅ'C���|�`  C1ׂ&��M�|sΊ�  C1����[�|h�  C2�BrM	�|\O��  C253�ޓ��|P��  C2X#�+�|DЃ   C2z���M�|9   C2����t��|-Q@  C2��	��2�|!��`  C2Ҕa���|�{�  C2��~���|
��  C3�5#��{�Sw�  C3�?%�{���  C3+��U�{��t   C3A��R��{��   C3Oy�?���{�Up@  C3Wy��h�{Õ�`  C3bF��4�{��l�  C3k�����{��  C3t�M�{�Wh�  C3~Ő}�{����  C3���W>��{��e   C3��
�ܦ�{}�   C3�b��Œ�{qYa@  C3�k%ң��{e��`  C3��K칄�{Y�]�  C4i�+��{N۠  C4"p!v7�{B[Y�  C42�R`z��{6���  C4>��|�e�{*�V   C4I B�{�   C4OM(q��{]R@  C4V��'to�{��`  C4g���8��z��N�  C4����])�z�̠  C4��01�z�_J�  C4�!�I4J�z؟��  C4�:�v���z��G   C4ܸ�H���z� �   C4��>&��z�aC@  C3kŐS>��z���`  C2���>��z��?�  C1�x��`��z�"��  C1v7!���z�c;�  C1]��c��zz���  C0��6�W��zn�8   C0��x�8�zc$�   C0� ���zWe4@  C0�=�;��zK��`  C0rus��z?�0�  C0bR�f*��z4&��  C0U�Jwq�z(g,�  C0[�~f��z���  C0t6:ϖ�z�)   C0�jYM�z(�   C0�q���&�y�i%@  C0�@*�y���`  C1�̟���y��!�  C1.�%����y�*��  C1CX��'�y�k�  C1Xѣ���y����  C1j�j��y��   C1xj*�AW�y�,�   C1��Z��U�y�m@  C1���՟�y���`  C1é����y���  C1�%
���yx.��  C2�;�1M�ylo�  C2&�_��y`���  C2A��W�yT�   C2W����yI0�   C2k�lH���y=q@  C2{�=	F)�y1��`  C2�`'�f�y%��  C2�o[o	�y2��  C2�#�׶�yr��  C2�������y�}�  C2�٣?��x���   C2�����x�4z   C3�{#���x�t�@  C35� `��xӵv`  C3P�)�,��x���  C3i��S��x�6r�  C3�Z�x�v��  C3��-0j��x��n�  C3�dɕq��x���   C3��2���x�8k   C3�_!̐>�x�x�@  C3��st
�xu�g`  C3��S�xi��  C3��8A�x^:c�  C3��<Y��xRz��  C3�3�"���xF�_�  C4 ۊ�X�x:��   C4�rW{�x/<\   C4;۟f�x#|�@  C4S������x�X`  C4ih�}��x�ր  C4zx�q���x >T�  C4�C["���w�~��  C4�h'w���w�P�  C4������w���   C4���>��w�@M   C4�l'���wŀ�@  C4�T<*��w��I`  C4ܨ}����w�ǀ  C3�)R0��w�BE�  C2�c��w��w����  C1���Kz�w��A�  C1g�彾�w�   C12E��wsD>   C0��}���wg��@  C0�*��2J�w[�:`  C0tv+�,�wP��  C0S��2��wDF6�  C0;\���v�w8���  C02=�z�w,�2�  C04�_0���w!�   C0FA))L��wH/   C0hN�͵��w	��@  C0��B���v��+`  C0�5.+��v�	��  C0ٸS�k�v�J'�  C0�V�,��vڊ��  C1������v��#�  C12�#��y�v��   C1I��)���v�L    C1Z\t�.o�v���@  C1l�r��v��`  C1~�P٭��v���  C1��D�<�v�N�  C1����̕�v|���  C1���ŗ;�vp��  C1�#�$&��ve�   C1�W���@�vYP   C1�m���vM��@  C1���T��vA�`  C2!Үj��v6��  C2?L����v*R	�  C2Yl��`;�v���  C2o��M�v��  C2�t����v�   C2��q�<;�u�T   C2�Kg��;�u@  C2��G*��u���`  C2�����u�|�  C2�e�)��u�U��  C2�جE[��u��x�  C31#s*��u����  C3:��O>�u�u   C3U�s�!��u�W�   C3n`n �"�u��q@  C3�`�5�V�u���`  C3�^e? ��uzm�  C3��˪�N�unY�  C3��B^y8�ub�i�  C3�ν.�uV���  C3����K��uKf   C3ӛ�}���u?[�   C3ۖ �YK�u3�b@  C3�C�w�u'��`  C3�L:%��u^�  C4c� �e�u]ܠ  C44f~۝l�u�Z�  C4O�Zy0�t����  C4g[�����t�W   C4|*5�t�_�   C4��ӕ��tՠS@  C4Zy��� �t���`  C31I=����t�!O�  C2l��m�A�t�a͠  C1�C)*�k�t��K�  C16PW�#��t����  C0нӶS�t�#H   C0�>����t�c�   C0g�]*p�tw�D@  C0M���m�tk��`  C0:q3&�t`%@�  C0.���1�tTe��  C0$pUc7��tH�<�  C00A}���t<��  C0�o��v�t1'9   C0�<�'b�t%g�   C0��y6��t�5@  C03��fG_�t�`  C0J�G�{!�t)1�  C0Y�K�ܧ�s�i��  C0h�����s�-�  C0���g-�s���  C0���[���s�+*   C0� �:z�s�k�   C0�(�?�s��&@  C1@ｱ��s��`  C11w���s�-"�  C1L���P��s�m��  C1f5-��}�s���  C1|�Ck�l�s���  C1����f�su/   C1�� C��sio�   C1����s]�@  C1��);~��sQ�`  C1�U$��sF1�  C1٧Ȧ��s:q��  C1����� �s.��  C26F�lY�s"��  C27ţ|�s3   C2Y�p�F�ss�   C2x)��E�r��@  C2�����r��`  C2�qL[7�r�5�  C2ǈ(j�s�r�u��  C2�*H ��rж �  C2��h�j��r��~�  C2���X��r�6�   C3ng����r�w{   C3��%<�r���@  C3%�o����r��w`  C3/g�
+��r�8��  C3:�n޺�r~ys�  C3P2Ɨ�n�rr���  C3p�����rf�o�  C3�
�sA�r[:�   C3����b�rO{l   C3�^@�rC��@  C3�ATQ��r7�h`  C3�Dll���r,<�  C4��Ψ�r }d�  C4�D��]�r���  C4(����_�r�`�  C43$�#,��q�>�   C4;V�.�q�]   C4A�B�4�q��@  C4G��q� Y`  C4M��B�q�@׀  C4_Ɩ���qU�  C4| �M0_�q����  C4�4꼱R�q�Q�  C4~P�3��q�B�   C3>��.U`�q��N   C2\��G�q���@  C1��j��@�q|J`  C1�k��!�qpDȀ  C0���- �qd�F�  C0_���qX���  C0(U�\mh�qMB�  C/�1?z8�qAF�   C/��рF��q5�?   C/������q)ǽ@  C/��2X[m�q;`  C/x�o����qH��  C/t�ed��q�7�  C/_�^���p�ɵ�  C/O-�Z���p�
3�  C/5�ωt�p�J�   C/}H�0[�p׋0   C/?HjEp�p�ˮ@  C/p��G��p�,`  C/�es1���p�L��  C0��d9��p��(�  C0+Sa�~��p�ͦ�  C0R��*�p�$�  C0v��MGx�p�N�   C0�S��bM�py�!   C0������pmϟ@  C0��'�	�pb`  C0��8��T�pVP��  C1m)�Ka�pJ��  C1���p>ї�  C1(Jn�)P�p3�  C18��C�p'R�   C1F�+{���p�   C1T�A�h��pӐ@  C1c^���\�p`  C1|�Z86�o�   C1�2X�o�*@  C1ğ��6}�o���  C1�*p'	��o�,�  C2~&5{��o��
   C2%�2P�(�o{.@  C2@�1ڐ�oc��  C2Y�����oL/��  C2o��@���o4��   C2�W����o1�@  C2��0]�W�o��  C2������n�3��  C2�v K�nִ�   C2���ك*�n�5�@  C2��%��x�n���  C2ьq�ԃ�n�7��  C2��J���nx��   C3	)�4*��na9�@  C3)�� �nI�Հ  C3G������n2;��  C3b�M���n��   C3x�-�k�n=�@  C3��w�~S�m�ƀ  C3�Qm�)�m�?��  C3�|#*�g�m���   C3�!�S8�m�A�@  C3������m�·�  C3�dz��!�mvC��  C3��:l��m^İ   C3������mGE�@  C4{�G���m/ƨ�  C4)?����mG��  C4<�q��O�m ȡ   C4Kl'�l�I�@  C3�1�n{�l�ʙ�  C2����Z �l�K��  C2'�P��l�̒   C1q�n{���l�M�@  C0�.�B�\�lsΊ�  C0��r�o��l\O��  C0\KA����lDЃ   C0.����B�l-Q@  C0�a�|��l�{�  C/�X|6��k�Sw�  C/�?4fS�k��t   C/�y����k�Up@  C/�0TU��k��l�  C/���w+0�k�Wh�  C/p�ѯV��k��e   C/O$j�O�kqYa@  C/,�����kY�]�  C/S����kB[Y�  C.�0"��c�k*�V   C/b9�}��k]R@  C/8�f� ��j��N�  C/h �E��j�_J�  C/� v��L�j��G   C/�S�< X�j�aC@  C/������j��?�  C0 P����j�c;�  C0&��jn�8   C0J�3���jWe4@  C0k��L��j?�0�  C0�DEhD[�j(g,�  C0�⋛���j�)   C0�҃30a�i�i%@  C0ǀT���i��!�  C0�0�)��i�k�  C0���c);�i��   C0�y"��}�i�m@  C1'{�&�i���  C1F�	���ilo�  C1g�s�g�iT�   C1���&M�i=q@  C1����>S�i%��  C1�U����ir��  C1ؖ�HL��h���   C1�E��p�h�t�@  C2��%���h���  C2^�b�b�h�v��  C2 �F�q��h���   C2-dGT���h�x�@  C29	���L�hi��  C2E�29L��hRz��  C2]�̓��h:��   C2�.��1�h#|�@  C2�0�Q��h�ր  C2�V�Wt��g�~��  C2��蛙�g���   C2����G�gŀ�@  C3M�#�f�g�ǀ  C3 �����g����  C3!9��7�g�   C3+����p�gg��@  C35j~ұ�gP��  C3J^�[X��g8���  C3iW���g!�   C3�j�d�g	��@  C3���W�)�f�	��  C3Ƹ�4"�fڊ��  C3�ˎsW}�f��   C4����f���@  C4C �,��f���  C4�l���f|���  C3)���,�fe�   C2_e^\�z�fM��@  C1�� ����f6��  C1���{�f���  C0�Q�,�f�   C0JSMf���e@  C0���h��e�|�  C/��T��e��x�  C/e���g�e�u   C/$M]���e��q@  C.�A�vx�ezm�  C.�B1�|�eb�i�  C.�Z���~�eKf   C.�~��B�e3�b@  C/O���e^�  C/2=�����e�Z�  C/3����d�W   C/1t�}�dՠS@  C/*�ה���d�!O�  C/$�"��o�d��K�  C/5_#;��d�#H   C/d��p��dw�D@  C/�9=j�d`%@�  C/���^�|�dH�<�  C/�Z�ů�d1'9   C0�ֈ<�d�5@  C0>�h���d)1�  C0b�r����c�-�  C0����A��c�+*   C0���.'�c��&@  C0ʊ
\���c�-"�  C0޶ uu�c���  C0��Y��o�cu/   C1#�6J��c]�@  C1�s���cF1�  C1,�%���c.��  C1;.�M��c3   C1H� W�b��@  C1Ve�S}�b�5�  C1o�l����bж �  C1�k����b�6�   C1�F}�}��b���@  C1ل��ۓ�b�8��  C1��Gg��br���  C2�Y2��b[:�   C27�v/��bC��@  C2Rm�ݟ�b,<�  C2kU.(@�b���  C2�'lpy�a�>�   C2��F���a��@  C2��p��!�a�@׀  C2�9���a����  C2��
M4�a�B�   C2թS�y�a���@  C2����U�apDȀ  C2�X����aX���  C2�B���
�aAF�   C3/S�N�a)ǽ@  C3�:�aH��  C35�#��@�`�ɵ�  C3U�f���`�J�   C3u�����`�ˮ@  C3�����R�`�L��  C3�g��%�`�ͦ�  C3�v��u�`�N�   C3���Z �`mϟ@  C3�{�;�`VP��  C4�����`>ї�  C4�I��`'R�   C4%8԰��`Ӑ@  C3��\0R��_�   C2�Ő!W�_���  C1�P�E4��_��
   C1<��6n��_c��  C0��y���_4��   C0T�0���_��  C0׊f�G�^ִ�   C/ʰ(��&�^���  C/�QQ� Y�^x��   C/f|��d��^I�Հ  C/I^�ʁ�^��   C/)<���T�]�ƀ  C/�hx��]���   C/���A�]�·�  C/��4S��]^İ   C.�|�c�]/ƨ�  C.���6���] ȡ   C.�����9�\�ʙ�  C.��E?��\�̒   C.���c��\sΊ�  C.���5�}�\DЃ   C/ ��複�\�{�  C/U[�x�[��t   C/�#��ac�[��l�  C/��!Hc�[��e   C0 ��a���[Y�]�  C0=�:F�t�[*�V   C0]43q�B�Z��N�  C0{#K;Oo�Z��G   C0�����L�Z��?�  C0��D�s�Zn�8   C0������Z?�0�  C0�����Z�)   C0��6$3�Y��!�  C0���#eH�Y��   C1W}t�Y���  C16D�ރ�YT�   C1'�F���Y%��  C16S1*��X���   C1O�w+��X���  C1s7
��X���   C1�������Xi��  C1�5�X�X:��   C1֗�����X�ր  C1��1���W���   C2�~�kq�W�ǀ  C2(<�A���W�   C22і~aV�WP��  C2A\����W!�   C2NC�E~��V�	��  C2Y{�qD�V��   C2dc��W��V���  C2{6��,�Ve�   C2�D����V6��  C2��߄�!�V�   C2�Lp���U�|�  C2��4�|��U�u   C3=�1͉�Uzm�  C3(���v*�UKf   C3;���;��U^�  C3K�r���T�W   C3YX�/��T�!O�  C3c�?g���T�#H   C3mf��o��T`%@�  C3vz�d��T1'9   C3����rw�T)1�  C3�e� X.�S�+*   C3�}b�	�S�-"�  C3�Ƹ���Su/   C4���/��SF1�  C4�bv�y�S3   C40����R�5�  C3������R�6�   C2�IX�A�R�8��  C2	(cWT��R[:�   C1M7���R,<�  C0���;�Q�Q�>�   C0\�(,y�Q�@׀  C0`�����Q�B�   C/ÓBQG��QpDȀ  C/p��s���QAF�   C/(R����QH��  C.�d+��T�P�J�   C.���3V�P�L��  C.�h�n3��P�N�   C.�����PVP��  C.��htSg�P'R�   C.���ðY�O�   C.����4��O��
   C.�g�lo��O4��   C.��v�`Z�Nִ�   C.�l�g���Nx��   C.��*H��N��   C/&*R����M���   C/N������M^İ   C/kVFp�M ȡ   C/�M�`�w�L�̒   C/�%��7~�LDЃ   C0"�����K��t   C0L�T�.��K��e   C0m�6F*�K*�V   C0�`F�ig�J��G   C0��]���Jn�8   C0�Df��J�)   C0⓿'�P�I��   C0��.�&`�IT�   C1��M�H���   C19��8�H���   C10���h��H:��   C1?�A8���G���   C1M��:�G�   C1Y�G���G!�   C1ht�h�F��   C1�bx����Fe�   C1���t��F�   C1�R1.�"�E�u   C1��I�eG�EKf   C2�k�\��D�W   C2++6��D�#H   C2KD��`/�D1'9   C2eam�k��C�+*   C2v���Sj�Cu/   C2�8����C3   C2��`Y��B�6�   C2��iB;�B[:�   C2�ɸ"x�A�>�   C2��l�T�A�B�   C2��۷��AAF�   C2ܕ�#�%�@�J�   C2�U��d��@�N�   C2������@'R�   C3�Y��?��
   C3>ՠ� �>ִ�   C3\�z7;��>��   C3w�aۏ,�=^İ   C3�L䀀a�<�̒   C3��ľ}��;��t   C3�G�y��;*�V   C3�!$W��:n�8   C3�PR�(f�9��   C3���y�8���   C3��p2y�8:��   C3���7�   C3���y"J�6��   C4:{b
��6�   C4&����5Kf   C4B��xy��4�#H   C4Sٓ1*�3�+*   C3m�t�0��33   C2s�����2[:�   C1�>��0�1�B�   C0�<R��0�J�   C0�����0'R�   C09���/<�.ִ�   C/�Fmu��-^İ   C/�VC�{��+��t   C/X�,y��*n�8   C/"�P��(���   C.��Nt��'�   C.�5�g�Q�&�   C.� �"�$�#H   C.��>�C�#3   C.ⴙe+y�!�B�   C/ ݫ�>l� 'R�   C/�r�^İ   C/84�����n�8   C/M5k�x��   C/f	4�*���#H   C/�m��2���B�   C/��}D�^İ   C05��ȓ��   C0)��A��B�   C0B�ޫ���   C0V�}�����   C0i�l�        C0{�H��IA��   C0�;@�A��   C0�"�P\QB�B�   C0�Tr�(�B�   C0�Dt��B^İ   C0�+0lwqB�B�   C0�f�&~VB�#H   C17&�[�B�   C1;�ޕ7\Bn�8   C1\��D��B^İ   C1{�O�oB 'R�   C1���
QB!�B�   C1�̟�UEB#3   C1�|M\_�B$�#H   C1�-Gс�B&�   C1�����B'�   C2
K5 YB(���   C2�!/�B*n�8   C2(�Cw��B+��t   C26"��B-^İ   C2A�=h�5B.ִ�   C2L�j�y�B0'R�   C2Y�p�ӞB0�J�   C2p��]5B1�B�   C2���ʰB2[:�   C2�0|O�B33   C2�v����B3�+*   C2������B4�#H   C3	��RBJB5Kf   C3v`\C�B6�   C31�@��B6��   C3B^AC]�B7�   C3O!��"B8:��   C3Z. DZ�B8���   C3cI:SB9��   C3l�5�YrB:n�8   C3�zva�OB;*�V   C3������B;��t   C3ʳwU��B<�̒   C3���XB=^İ   C4 A�'��B>��   C4�թ��B>ִ�   C4&�I
�SB?��
   C4-�-��B@'R�   C3VUΈB@�N�   C2k���B@�J�   C1��Y�~�BAAF�   C1�X�n�BA�B�   C0�G��&�BA�>�   C03��;�BB[:�   C/�_:zBB�6�   C/�*��BC3   C/q�9�&tBCu/   C/?��lRBC�+*   C/#�ΥBD1'9   C/�z~vBD�#H   C/
F���BD�W   C.�T�ރuBEKf   C.�V��RfBE�u   C.�=�H�0BF�   C.�����BFe�   C.�%���BF��   C.{�C��BG!�   C.���ZABG�   C.�ҐX�!BG���   C/���^BH:��   C/* �j^<BH���   C/L��=1�BH���   C/rc*J�BIT�   C/��nڞBI��   C/�MݭݟBJ�)   C0$!��![BJn�8   C0F�@q��BJ��G   C0j��w<%BK*�V   C0�s�h[BK��e   C0��I��}BK��t   C0��R�zRBLDЃ   C0ްt��pBL�̒   C0�m��N�BM ȡ   C0�}BV31BM^İ   C1CD�:FBM���   C1X��аBN��   C1+rQ�BNx��   C1H	��b�BNִ�   C1c��ԀWBO4��   C1����#!BO��
   C1��ϑ%BO�   C1�>֕�BP'R�   C1�s����BPVP��  C2
�^���BP�N�   C2&�Y���BP�L��  C2A#=�y�BP�J�   C2Y'G�EBQH��  C2ok`K�BQAF�   C2�I^�]BQpDȀ  C2�f!�9�BQ�B�   C2���[I�BQ�@׀  C2�z��ЅBQ�>�   C2�CB�lBR,<�  C2̈�x�	BR[:�   C2���.�BR�8��  C2���b�BR�6�   C2�^t��BR�5�  C3%�O��BS3   C3"+5C"�BSF1�  C3B�h���BSu/   C3a�*fݣBS�-"�  C3~����KBS�+*   C3��5c�BT)1�  C3���ջOBT1'9   C3��{��BT`%@�  C3�e�92�BT�#H   C3���BT�!O�  C4�&_gBT�W   C4Ew��iBU^�  C4�z4�:BUKf   C3~�{3�BUzm�  C2����	BU�u   C1���`�-BU�|�  C10z2 �BV�   C0��x%��BV6��  C03�F��LBVe�   C/�����UBV���  C/���G�BV��   C/l2�X�BV�	��  C/;Mx��BW!�   C/��9��BWP��  C.�i�H�BW�   C.����dIBW�ǀ  C.�2�2/�BW���   C.��X-�BX�ր  C.u���BX:��   C.g�
S�BXi��  C.n���BX���   C.w�1�BX���  C.x����XBX���   C.{�Z!BY%��  C.������BYT�   C.�zd�N�BY���  C/)��c��BY��   C/n��Zh�BY��!�  C/�<�`BZ�)   C/ػS)�BZ?�0�  C0����BZn�8   C0!tj%5BZ��?�  C0)����BZ��G   C09s���BZ��N�  C0HڨCWB[*�V   C0\
��)B[Y�]�  C0zM�,v�B[��e   C0��t�pB[��l�  C0��?ðGB[��t   C0�dw��dB\�{�  C1����B\DЃ   C1'�U#8�B\sΊ�  C1Gb��l�B\�̒   C1c«
��B\�ʙ�  C1}��4��B] ȡ   C1��P2�B]/ƨ�  C1���.�B]^İ   C1�Qk���B]�·�  C1�@����B]���   C1�y�C�B]�ƀ  C1�&�Έ�B^��   C2G��»B^I�Հ  C2��T�5B^x��   C2�M��XB^���  C2',�|�B^ִ�   C24��bhB_��  C2L��t:�B_4��   C2mƂ�y�B_c��  C2��w�rB_��
   C2���?�WB_���  C2�s3��)B_�   C2�B�_��B`Ӑ@  C3fn�DB`'R�   C3�{}�B`>ї�  C34Y�4a�B`VP��  C3Mӕ3B`mϟ@  C3b���B`�N�   C3g)�,�B`�ͦ�  C3u�}ՔB`�L��  C3�+�^1�B`�ˮ@  C3����NB`�J�   C3�*J�B`�ɵ�  C3�Q�?��BaH��  C3��w�Ba)ǽ@  C3���o��BaAF�   C3��t�BaX���  C4
j�wjIBapDȀ  C4Є��]Ba���@  C3;��tQBa�B�   C2[B�W��Ba����  C1������Ba�@׀  C0��M���Ba��@  C0nN�ؘBa�>�   C0�b�Bb���  C/��wF�Bb,<�  C/���!7.BbC��@  C/C��ȋ4Bb[:�   C/�ӡ{�Bbr���  C.���O@Bb�8��  C.�R�z(CBb���@  C.ǭ蘢�Bb�6�   C.��w���Bbж �  C.��"�`Bb�5�  C.t��}L+Bb��@  C.uW�.nwBc3   C.�M� Bc.��  C.�ѷD��BcF1�  C.ϗ���Bc]�@  C/�l�S�Bcu/   C/a�lQR!Bc���  C/����Bc�-"�  C/؃��$�Bc��&@  C0��r1�Bc�+*   C0%/�,Bc�-�  C01�5�Bd)1�  C0C��guBd�5@  C0S����|Bd1'9   C0bL����BdH�<�  C0r�eBd`%@�  C0���_΃Bdw�D@  C0������Bd�#H   C0֖��E>Bd��K�  C0��zBd�!O�  C1��)�BdՠS@  C1>�}�Bd�W   C1]8�b��Be�Z�  C1y���Be^�  C1��¤~(Be3�b@  C1��vA�BeKf   C1���'��Beb�i�  C1�>���Bezm�  C1���|��Be��q@  C1���8�Be�u   C2�o[�Be��x�  C2��