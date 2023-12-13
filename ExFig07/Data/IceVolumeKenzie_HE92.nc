CDF  �   
      time       lon       lat          	   CDI       @Climate Data Interface version 2.0.3 (https://mpimet.mpg.de/cdi)   Conventions       CF-1.5     source        #PISM development no-version-control    command      / /home/m/m300792/ModelCodes/mPISM_build/install_03-10-22/bin/pismr -i pism_-065000/pism_-065000.nc -stress_balance ssa+sia -y 100 -ys -65000 -o pism_-064900/pism_-064900.nc -surface_given_file forcing/smb_-064900.nc -surface given -ocean_th_file forcing/oce_-064900.nc -ocean th -topg_file pism_-065000/pism_-065000_new_topg.nc -bed_def prescribed -uplift_file pism_-065000/pism_-065000_dbdt.nc -extra_times 10 -extra_vars thk,usurf,bheatflx,velsurf_mag,velbase_mag,tillwat,bmelt,velbase,flux_divergence,velsurf,velbar,wvelsurf,wvelbase,climatic_mass_balance_cumulative,potential_climatic_mass_balance_cumulative,tempbase,tempsurf,mask,temppabase,tempicethk_basal,topg,nonneg_flux_cumulative,grounded_basal_flux_cumulative,floating_basal_flux_cumulative,discharge_flux_cumulative,Href,taud_mag,lat,lon,calving_mask,sea_level -extra_file pism_-064900/pism_-064900_extra -extra_split -ts_times yearly -ts_file pism_-064900/pism_-064900_ts.nc -topg_to_phi 15.0,30.0,-300.0,400.0 -ssafd_ksp_rtol 1e-7 -enthalpy_temperate_conductivity_ratio 0.001 -pseudo_plastic -pseudo_plastic_q .25 -pseudo_plastic_uthreshold 70 -till_effective_fraction_overburden 0.04 -tauc_slippery_grounding_lines -sia_e 8 -ssa_e 1 -ssa_upwind_factor 0.7 -calving eigen_calving,thickness_calving -thickness_calving_threshold 100 -part_redist -part_grid -cfbc -kill_icebergs -eigen_calving_K 1e16 -verbose 2 -sliding_scale_factor_reduces_tauc 1 -th_heat_coefficient 4e-5 -th_salinity_coefficient 4e-8 -o_format netcdf3 -sediment_file till_mask.nc -no_divide_tillphi -periodicity none -ocean_th_period 1 -o_order yxz
    proj4         j+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs     history_of_appended_files        �Fri Jul 29 07:49:41 2022: Appended file outFinal2 had following "history" attribute:
Fri Jul 29 07:49:39 2022: ncap2 -O -s lon_bnds=float(lon_bnds) outFinal1 outFinal2
Fri Jul 29 07:49:39 2022: ncap2 -O -s lat_bnds=float(lat_bnds) outFinal outFinal1
Fri Jul 29 07:49:38 2022: cdo -add thkNew thkOld outFinal
Fri Jul 29 07:49:36 2022: cdo -mul -selvar,thk pism_-064900/pism_-064900.nc MaskLaurentide thkNew
     NCO       _netCDF Operators version 5.0.6 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)    history      7�Wed Sep 14 12:02:51 2022: cdo -mergetime File000001.nc File000002.nc File000003.nc File000004.nc File000005.nc File000006.nc File000007.nc File000008.nc File000009.nc File000010.nc File000011.nc File000012.nc File000013.nc File000014.nc File000015.nc File000016.nc File000017.nc File000018.nc File000019.nc File000020.nc File000021.nc File000022.nc File000023.nc File000024.nc File000025.nc File000026.nc File000027.nc File000028.nc File000029.nc File000030.nc File000031.nc File000032.nc File000033.nc File000034.nc File000035.nc File000036.nc File000037.nc File000038.nc File000039.nc File000040.nc File000041.nc File000042.nc File000043.nc File000044.nc File000045.nc File000046.nc File000047.nc File000048.nc File000049.nc File000050.nc File000051.nc File000052.nc File000053.nc File000054.nc File000055.nc File000056.nc File000057.nc File000058.nc File000059.nc File000060.nc File000061.nc File000062.nc File000063.nc File000064.nc File000065.nc File000066.nc File000067.nc File000068.nc File000069.nc File000070.nc File000071.nc File000072.nc File000073.nc File000074.nc File000075.nc File000076.nc File000077.nc File000078.nc File000079.nc File000080.nc File000081.nc File000082.nc File000083.nc File000084.nc File000085.nc File000086.nc File000087.nc File000088.nc File000089.nc File000090.nc File000091.nc File000092.nc File000093.nc File000094.nc File000095.nc File000096.nc File000097.nc File000098.nc File000099.nc File000100.nc File000101.nc File000102.nc File000103.nc File000104.nc File000105.nc File000106.nc File000107.nc File000108.nc File000109.nc File000110.nc File000111.nc File000112.nc File000113.nc File000114.nc File000115.nc File000116.nc File000117.nc File000118.nc File000119.nc File000120.nc File000121.nc File000122.nc File000123.nc File000124.nc File000125.nc File000126.nc File000127.nc File000128.nc File000129.nc File000130.nc File000131.nc File000132.nc File000133.nc File000134.nc File000135.nc File000136.nc File000137.nc File000138.nc File000139.nc File000140.nc File000141.nc File000142.nc File000143.nc File000144.nc File000145.nc File000146.nc File000147.nc File000148.nc File000149.nc File000150.nc File000151.nc File000152.nc File000153.nc File000154.nc File000155.nc File000156.nc File000157.nc File000158.nc File000159.nc File000160.nc File000161.nc File000162.nc File000163.nc File000164.nc File000165.nc File000166.nc File000167.nc File000168.nc File000169.nc File000170.nc File000171.nc File000172.nc File000173.nc File000174.nc File000175.nc File000176.nc File000177.nc File000178.nc File000179.nc File000180.nc File000181.nc File000182.nc File000183.nc File000184.nc File000185.nc File000186.nc File000187.nc File000188.nc File000189.nc File000190.nc File000191.nc File000192.nc File000193.nc File000194.nc File000195.nc File000196.nc File000197.nc File000198.nc File000199.nc File000200.nc File000201.nc File000202.nc File000203.nc File000204.nc File000205.nc File000206.nc File000207.nc File000208.nc File000209.nc File000210.nc File000211.nc File000212.nc File000213.nc File000214.nc File000215.nc File000216.nc File000217.nc File000218.nc File000219.nc File000220.nc File000221.nc File000222.nc File000223.nc File000224.nc File000225.nc File000226.nc File000227.nc File000228.nc File000229.nc File000230.nc File000231.nc File000232.nc File000233.nc File000234.nc File000235.nc File000236.nc File000237.nc File000238.nc File000239.nc File000240.nc File000241.nc File000242.nc File000243.nc File000244.nc File000245.nc File000246.nc File000247.nc File000248.nc File000249.nc File000250.nc File000251.nc File000252.nc File000253.nc File000254.nc File000255.nc File000256.nc File000257.nc File000258.nc File000259.nc File000260.nc File000261.nc File000262.nc File000263.nc File000264.nc File000265.nc File000266.nc File000267.nc File000268.nc File000269.nc File000270.nc File000271.nc File000272.nc File000273.nc File000274.nc File000275.nc File000276.nc File000277.nc File000278.nc File000279.nc File000280.nc File000281.nc File000282.nc File000283.nc File000284.nc File000285.nc File000286.nc File000287.nc File000288.nc File000289.nc File000290.nc File000291.nc File000292.nc File000293.nc File000294.nc File000295.nc File000296.nc File000297.nc File000298.nc File000299.nc File000300.nc File000301.nc File000302.nc File000303.nc File000304.nc File000305.nc File000306.nc File000307.nc File000308.nc File000309.nc File000310.nc File000311.nc File000312.nc File000313.nc File000314.nc File000315.nc File000316.nc File000317.nc File000318.nc File000319.nc File000320.nc File000321.nc File000322.nc File000323.nc File000324.nc File000325.nc File000326.nc File000327.nc File000328.nc File000329.nc File000330.nc File000331.nc File000332.nc File000333.nc File000334.nc File000335.nc File000336.nc File000337.nc File000338.nc File000339.nc File000340.nc File000341.nc File000342.nc File000343.nc File000344.nc File000345.nc File000346.nc File000347.nc File000348.nc File000349.nc File000350.nc File000351.nc File000352.nc File000353.nc File000354.nc File000355.nc File000356.nc File000357.nc File000358.nc File000359.nc File000360.nc File000361.nc File000362.nc File000363.nc File000364.nc File000365.nc File000366.nc File000367.nc File000368.nc File000369.nc File000370.nc File000371.nc File000372.nc File000373.nc File000374.nc File000375.nc File000376.nc File000377.nc File000378.nc File000379.nc File000380.nc File000381.nc File000382.nc File000383.nc File000384.nc File000385.nc File000386.nc File000387.nc File000388.nc File000389.nc File000390.nc File000391.nc File000392.nc File000393.nc File000394.nc File000395.nc File000396.nc File000397.nc File000398.nc File000399.nc File000400.nc File000401.nc File000402.nc File000403.nc File000404.nc File000405.nc File000406.nc File000407.nc File000408.nc File000409.nc File000410.nc File000411.nc File000412.nc File000413.nc File000414.nc File000415.nc File000416.nc File000417.nc File000418.nc File000419.nc File000420.nc File000421.nc File000422.nc File000423.nc File000424.nc File000425.nc File000426.nc File000427.nc File000428.nc File000429.nc File000430.nc File000431.nc File000432.nc File000433.nc File000434.nc File000435.nc File000436.nc File000437.nc File000438.nc File000439.nc File000440.nc File000441.nc File000442.nc File000443.nc File000444.nc File000445.nc File000446.nc File000447.nc File000448.nc File000449.nc File000450.nc File000451.nc File000452.nc File000453.nc File000454.nc File000455.nc File000456.nc File000457.nc File000458.nc File000459.nc File000460.nc File000461.nc File000462.nc File000463.nc File000464.nc File000465.nc File000466.nc File000467.nc File000468.nc File000469.nc File000470.nc File000471.nc File000472.nc File000473.nc File000474.nc File000475.nc File000476.nc File000477.nc File000478.nc File000479.nc File000480.nc File000481.nc File000482.nc File000483.nc File000484.nc File000485.nc File000486.nc File000487.nc File000488.nc File000489.nc File000490.nc File000491.nc File000492.nc File000493.nc File000494.nc File000495.nc File000496.nc File000497.nc File000498.nc File000499.nc File000500.nc File000501.nc File000502.nc File000503.nc File000504.nc File000505.nc File000506.nc File000507.nc File000508.nc File000509.nc File000510.nc File000511.nc File000512.nc File000513.nc File000514.nc File000515.nc File000516.nc File000517.nc File000518.nc File000519.nc File000520.nc File000521.nc File000522.nc File000523.nc File000524.nc File000525.nc File000526.nc File000527.nc File000528.nc File000529.nc File000530.nc File000531.nc File000532.nc File000533.nc File000534.nc File000535.nc File000536.nc File000537.nc File000538.nc File000539.nc File000540.nc File000541.nc File000542.nc File000543.nc File000544.nc File000545.nc File000546.nc File000547.nc File000548.nc File000549.nc File000550.nc File000551.nc File000552.nc File000553.nc File000554.nc File000555.nc File000556.nc File000557.nc File000558.nc File000559.nc File000560.nc File000561.nc File000562.nc File000563.nc File000564.nc File000565.nc File000566.nc File000567.nc File000568.nc File000569.nc File000570.nc File000571.nc File000572.nc File000573.nc File000574.nc File000575.nc File000576.nc File000577.nc File000578.nc File000579.nc File000580.nc File000581.nc File000582.nc File000583.nc File000584.nc File000585.nc File000586.nc File000587.nc File000588.nc File000589.nc File000590.nc File000591.nc File000592.nc File000593.nc File000594.nc File000595.nc File000596.nc File000597.nc File000598.nc File000599.nc File000600.nc File000601.nc File000602.nc File000603.nc File000604.nc File000605.nc File000606.nc File000607.nc File000608.nc File000609.nc File000610.nc File000611.nc File000612.nc File000613.nc File000614.nc File000615.nc File000616.nc File000617.nc File000618.nc File000619.nc File000620.nc File000621.nc File000622.nc File000623.nc File000624.nc File000625.nc File000626.nc File000627.nc File000628.nc File000629.nc File000630.nc File000631.nc File000632.nc File000633.nc File000634.nc File000635.nc File000636.nc File000637.nc File000638.nc File000639.nc File000640.nc File000641.nc File000642.nc File000643.nc File000644.nc File000645.nc File000646.nc File000647.nc File000648.nc File000649.nc File000650.nc File000651.nc File000652.nc File000653.nc File000654.nc File000655.nc File000656.nc File000657.nc File000658.nc File000659.nc File000660.nc File000661.nc File000662.nc File000663.nc File000664.nc File000665.nc File000666.nc File000667.nc File000668.nc File000669.nc File000670.nc File000671.nc File000672.nc File000673.nc File000674.nc File000675.nc File000676.nc File000677.nc File000678.nc File000679.nc File000680.nc File000681.nc File000682.nc File000683.nc File000684.nc File000685.nc File000686.nc File000687.nc File000688.nc File000689.nc File000690.nc File000691.nc File000692.nc File000693.nc File000694.nc File000695.nc File000696.nc File000697.nc File000698.nc File000699.nc File000700.nc File000701.nc File000702.nc File000703.nc File000704.nc File000705.nc File000706.nc File000707.nc File000708.nc File000709.nc File000710.nc File000711.nc File000712.nc File000713.nc File000714.nc File000715.nc File000716.nc File000717.nc File000718.nc File000719.nc File000720.nc File000721.nc File000722.nc File000723.nc File000724.nc File000725.nc File000726.nc File000727.nc File000728.nc File000729.nc File000730.nc File000731.nc File000732.nc File000733.nc File000734.nc File000735.nc File000736.nc File000737.nc File000738.nc File000739.nc File000740.nc File000741.nc File000742.nc File000743.nc File000744.nc File000745.nc File000746.nc File000747.nc File000748.nc File000749.nc File000750.nc File000751.nc File000752.nc File000753.nc File000754.nc File000755.nc File000756.nc File000757.nc File000758.nc File000759.nc File000760.nc File000761.nc File000762.nc File000763.nc File000764.nc File000765.nc File000766.nc File000767.nc File000768.nc File000769.nc File000770.nc File000771.nc File000772.nc File000773.nc File000774.nc File000775.nc File000776.nc File000777.nc File000778.nc File000779.nc File000780.nc File000781.nc File000782.nc File000783.nc File000784.nc File000785.nc File000786.nc File000787.nc File000788.nc File000789.nc File000790.nc File000791.nc File000792.nc File000793.nc File000794.nc File000795.nc File000796.nc File000797.nc File000798.nc File000799.nc File000800.nc File000801.nc File000802.nc File000803.nc File000804.nc File000805.nc File000806.nc File000807.nc File000808.nc File000809.nc File000810.nc File000811.nc File000812.nc File000813.nc File000814.nc File000815.nc File000816.nc File000817.nc File000818.nc File000819.nc File000820.nc File000821.nc File000822.nc File000823.nc File000824.nc File000825.nc File000826.nc File000827.nc File000828.nc File000829.nc File000830.nc File000831.nc File000832.nc File000833.nc File000834.nc File000835.nc File000836.nc File000837.nc File000838.nc File000839.nc File000840.nc File000841.nc File000842.nc File000843.nc File000844.nc File000845.nc File000846.nc File000847.nc File000848.nc File000849.nc File000850.nc File000851.nc File000852.nc File000853.nc File000854.nc File000855.nc File000856.nc File000857.nc File000858.nc File000859.nc File000860.nc File000861.nc File000862.nc File000863.nc File000864.nc File000865.nc File000866.nc File000867.nc File000868.nc File000869.nc File000870.nc File000871.nc File000872.nc File000873.nc File000874.nc File000875.nc File000876.nc File000877.nc File000878.nc File000879.nc File000880.nc File000881.nc File000882.nc File000883.nc File000884.nc File000885.nc File000886.nc File000887.nc File000888.nc File000889.nc File000890.nc File000891.nc File000892.nc File000893.nc File000894.nc File000895.nc File000896.nc File000897.nc File000898.nc File000899.nc File000900.nc File000901.nc File000902.nc File000903.nc File000904.nc File000905.nc File000906.nc File000907.nc File000908.nc File000909.nc File000910.nc File000911.nc File000912.nc File000913.nc File000914.nc File000915.nc File000916.nc File000917.nc File000918.nc File000919.nc File000920.nc File000921.nc File000922.nc File000923.nc File000924.nc File000925.nc File000926.nc File000927.nc File000928.nc File000929.nc File000930.nc File000931.nc File000932.nc File000933.nc File000934.nc File000935.nc File000936.nc File000937.nc File000938.nc File000939.nc File000940.nc File000941.nc File000942.nc File000943.nc File000944.nc File000945.nc File000946.nc File000947.nc File000948.nc File000949.nc File000950.nc File000951.nc File000952.nc File000953.nc File000954.nc File000955.nc File000956.nc File000957.nc File000958.nc File000959.nc File000960.nc File000961.nc File000962.nc File000963.nc File000964.nc File000965.nc File000966.nc File000967.nc File000968.nc File000969.nc File000970.nc File000971.nc File000972.nc File000973.nc File000974.nc File000975.nc File000976.nc File000977.nc File000978.nc File000979.nc File000980.nc File000981.nc File000982.nc File000983.nc File000984.nc File000985.nc File000986.nc File000987.nc File000988.nc File000989.nc File000990.nc File000991.nc File000992.nc File000993.nc File000994.nc File000995.nc File000996.nc File000997.nc File000998.nc File000999.nc File001000.nc IceVolumeKenzie.nc
Wed Sep 14 11:57:03 2022: cdo -fldsum -mulc,10000 -mulc,10000 -selindexbox,211,427,600,415 -selvar,thk /work/ba0989/m300792/HE_Runs/ExperimentsComposite//HE92//pism_-064900/pism_-064900.nc Tmp/File000001.nc     CDO       @Climate Data Operators version 2.0.3 (https://mpimet.mpg.de/cdo)         time                standard_name         time   	long_name         time   units         seconds since 1-1-1    calendar      365_day    axis      T               D�   lon                standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X               D�   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y               D�   thk                       standard_name         land_ice_thickness     	long_name         land ice thickness     units         m      pism_intent       model_state             D�                �}Ȁ@�  C-V��z��}���   C-Cj���}�=   C-*�b.-��}�A�@  C-$8�&���}��9`  C-8���!��}�·�  C-d�t��j�}�5�  C-�<����}vC��  C-����
�}j�1�  C-��ȌÕ�}^İ   C.?�����}S.   C.��(J	�}GE�@  C.�܅���};�*`  C/y����C�}/ƨ�  C0
�꠨�}$&�  C0J{�L��}G��  C0I.=#�S�}�"�  C00�C|��} ȡ   C0C�,�#�|�	   C0Q���=��|�I�@  C05���M�|݊`  C0M�f��|�ʙ�  C0b0�E:�|��  C0x�^^W�|�K��  C0��ݕ���|���  C0��V��|�̒   C0�7
����|�   C1��P+��|�M�@  C1C��c(�|�`  C1���y��|sΊ�  C1����d��|h�  C1Ǟ��̪�|\O��  C1�]�A�6�|P��  C1ɴ&�7�|DЃ   C1���+N��|9   C1g�����|-Q@  C1[1Sr"�|!��`  C1Vp;L���|�{�  C1U���8�|
��  C1b�L �F�{�Sw�  C1}��W�A�{���  C1��-�S��{��t   C1�Ҥ��{��   C2
cN��}�{�Up@  C2N�%Cm��{Õ�`  C2����{��l�  C2��5;��{��  C2u�sP
��{�Wh�  C2 ;���{����  C0��m��{��e   C.%���.(�{}�   C-��]�I�{qYa@  C,~l��`��{e��`  C,9� <���{Y�]�  C,/�ё���{N۠  C,�H�4@��{B[Y�  C,�m̥�z�{6���  C-cv�Nڣ�{*�V   C-�h� �{�   C.������{]R@  C/��&�8�{��`  C/#���j�z��N�  C.�
��u�z�̠  C.�ߍj�b�z�_J�  C.�$��/�z؟��  C/ ��j��z��G   C.�W��z� �   C.u�M<ck�z�aC@  C.f}��R�z���`  C.mkU��z��?�  C.�Zֶ2�z�"��  C.�L!����z�c;�  C/6�Z����zz���  C/���Ph��zn�8   C0*G�zr��zc$�   C0hL�X��zWe4@  C0f�����zK��`  C0O��z���z?�0�  C0b?�X�7�z4&��  C0x���&�z(g,�  C0����z���  C0X�_�t`�z�)   C0h��K�z(�   C0g�:nK�y�i%@  C0�O�nL�y���`  C0��Q7sA�y��!�  C0أ$���y�*��  C1	֡Ut�y�k�  C1F�E��y����  C1���K��y��   C1�d�k�y�,�   C1�ﺀ��y�m@  C1��h���y���`  C1ǭ��q\�y���  C1�'s��yx.��  C1�/)��W�ylo�  C1�K��0-�y`���  C1}oH��E�yT�   C1x')}�}�yI0�   C1}1ûtU�y=q@  C1��3C�y1��`  C1���R ��y%��  C1��gf?�y2��  C0^���yr��  C0 3��{�y�}�  C/���l�{�x���   C/mK����x�4z   C.k��'��x�t�@  C. 	��?��xӵv`  C--,
R�x���  C,d�	�wE�x�6r�  C,�����x�v��  C+̮��r��x��n�  C+��s�}�x���   C+���W�;�x�8k   C+������x�x�@  C,B��*�q�xu�g`  C,����+�xi��  C-F���sS�x^:c�  C-��A�=�xRz��  C.i�l���xF�_�  C.t�.����x:��   C.=�mwJ�x/<\   C.?��$>�x#|�@  C.E7A��^�x�X`  C.OEu`�x�ր  C.qr
7k��x >T�  C.g�7v�w�~��  C. G�{��w�P�  C.����w���   C.)�g��'�w�@M   C.d�(#��wŀ�@  C.ȾxW�w��I`  C/P��W^��w�ǀ  C/��TW�w�BE�  C03�AF��w����  C00��}t�w��A�  C0 G�!�w�   C0*��v���wsD>   C0<bJM$_�wg��@  C0t�l��w[�:`  C0U�-��	�wP��  C0Oxz��K�wDF6�  C0[��Ö�w8���  C0q��sE�w,�2�  C0���`2 �w!�   C0����r�wH/   C0�ir�2�w	��@  C1%����2�v��+`  C1f�P�C��v�	��  C1�EO�>D�v�J'�  C1��j���vڊ��  C1�ytL` �v��#�  C1����[��v��   C1Ɩ8Zz��v�L    C1�G�[��v���@  C1�+U!�&�v��`  C1q$���[�v���  C1c}suW��v�N�  C1`y��g��v|���  C134�9�vp��  C0I��,�/�ve�   C/��	��vYP   C/�����vM��@  C.�Bm���vA�`  C.��ȥ���v6��  C-�I�ɒD�v*R	�  C-pu�*5��v���  C-7P��"�v��  C-б=Ѧ�v�   C-	���!`�u�T   C,n��i9.�u@  C,"��n�u���`  C,��T�u�|�  C,,��q�+�u�U��  C,f�dRX�u��x�  C,�#y�+ �u����  C-9���"�u�u   C-��O��u�W�   C.e� aM��u��q@  C.�$ܐ$��u���`  C.�s��Zc�uzm�  C.�/i���unY�  C.�Z�(�ub�i�  C.��L�^r�uV���  C.�"�!s�uKf   C/(B^R��u?[�   C/>�,�}�u3�b@  C.��:7	f�u'��`  C.��v���u^�  C.�asW.]�u]ܠ  C.�G�{ȭ�u�Z�  C/&��t_�t����  C/�,G��A�t�W   C0#�|��t�_�   C0a�WF6C�tՠS@  C0^��L��t���`  C0F��y"A�t�!O�  C0X�"l���t�a͠  C0i����T�t��K�  C0�\���5�t����  C0�(r��>�t�#H   C0h�ijP9�t�c�   C/�����tw�D@  C/4��gv��tk��`  C.���B���t`%@�  C.K7�5y<�tTe��  C.5f�8v�tH�<�  C.G"��~�t<��  C.m
����t1'9   C.��TaV\�t%g�   C.?��k���t�5@  C-�s��q�t�`  C-�����T�t)1�  C-_~�4�s�i��  C-7�o�/c�s�-�  C-*����Y�s���  C,��sL^��s�+*   C,C��o͒�s�k�   C,o�h���s��&@  C,p#��e�s��`  C,WF�a�s�-"�  C,��+�j�s�m��  C-ciǋ	��s���  C-����'��s���  C.��[2���su/   C.��@CQ�sio�   C.U�N\S�s]�@  C.W�x:��sQ�`  C.f׀�m��sF1�  C.v�/)���s:q��  C.��gû��s.��  C.ۭ����s"��  C.�������s3   C.��@�ss�   C.����,�r��@  C.��VУ�r��`  C/aM�ݚ�r�5�  C/���],�r�u��  C0E��7�rж �  C0^&΁��r��~�  C0Z����8�r�6�   C0C��jv�r�w{   C0V��X�r���@  C0k�J��r��w`  C0�kW�@8�r�8��  C/�΢���r~ys�  C.�z_��rr���  C.H��A��rf�o�  C-ܔ�?N�r[:�   C-�������rO{l   C-`�k�A�rC��@  C-Wb��t�r7�h`  C-o�Q����r,<�  C-����y��r }d�  C-���q
1�r���  C-�V�18�r�`�  C-\��cӜ�q�>�   C-E,t��q�]   C-:��O��q��@  C-8(�GD�q� Y`  C-M�ɧb�q�@׀  C-X�����qU�  C,�?�,.�q����  C,�l"Bl��q�Q�  C,���1���q�B�   C-/�vϡ��q��N   C-���o���q���@  C.���0�q|J`  C.�tP��C�qpDȀ  C/R�1��qd�F�  C.�qI� �qX���  C.���DM9�qMB�  C.��G%�w�qAF�   C.��B����q5�?   C.�wٻ���q)ǽ@  C.���z-.�q;`  C/Ot���qH��  C.�hݕ��q�7�  C.D���p�ɵ�  C.(�։��p�
3�  C.$6�b=��p�J�   C.Z+p����p׋0   C.��w����p�ˮ@  C.褍�@�p�,`  C/ �X���p�L��  C.�=L9���p��(�  C.���:q�p�ͦ�  C.S�m��X�p�$�  C.<BpD���p�N�   C.:�%��py�!   C-���
͠�pmϟ@  C->�~��]�pb`  C,�`*��9�pVP��  C,�]�傕�pJ��  C,͂�:��p>ї�  C,�P�t!�p3�  C-&�2e��p'R�   C-r���q��p�   C-��kQ��pӐ@  C.b�<t�p`  C-�3j0-��o�   C-y���pi�o�*@  C-\�ՙlZ�o���  C-S�l�l��o�,�  C-Iq��J\�o��
   C-Kl�=���o{.@  C,�`q�{�oc��  C,kAI�1�oL/��  C,T�J ��o4��   C,XEH��"�o1�@  C,�jq�o��  C,��i���n�3��  C-"�X[�nִ�   C-xsa�<��n�5�@  C-�aQY�j�n���  C-�/��#��n�7��  C-Oqz����nx��   C-6~��-��na9�@  C-070��nI�Հ  C-)�c����n2;��  C-;%;����n��   C-RE�A��n=�@  C,�˟o�m�ƀ  C,��u��:�m�?��  C,�\"��m���   C,�\(Ib��m�A�@  C-G8su�m�·�  C-��=���mvC��  C.#��1�m^İ   C.�`��1�mGE�@  C.q+�����m/ƨ�  C.5C�e�0�mG��  C.;G�����m ȡ   C.H#U/B�l�I�@  C.O�J��C�l�ʙ�  C.w?.���l�K��  C.�p��_*�l�̒   C.(8v�� �l�M�@  C-�+\�e�lsΊ�  C-���S��l\O��  C-�ov�m}�lDЃ   C-�"ʌ8��l-Q@  C./O?�$(�l�{�  C.f�_K�k�Sw�  C.�ص?���k��t   C.S�^�L��k�Up@  C-�R�`l��k��l�  C-���W�2�k�Wh�  C-�Z2����k��e   C-kez�^�kqYa@  C-Y�G.���kY�]�  C,�G����kB[Y�  C,�DP��S�k*�V   C,�IX���k]R@  C,���;��j��N�  C,���]�;�j�_J�  C,�\j���j��G   C-;�I��%�j�aC@  C-�2�i���j��?�  C-�����j�c;�  C-���.i�jn�8   C-a��j �jWe4@  C-L2���j?�0�  C-I�����j(g,�  C-I������j�)   C-I
h��#�i�i%@  C-$���i��!�  C,�ȓ�)g�i�k�  C,d�?~��i��   C,_��i�{�i�m@  C,�����i���  C,���d�ilo�  C-4"��R��iT�   C-�#���i=q@  C-���~���i%��  C-����W�ir��  C-B�b�� �h���   C- 4�.�h�t�@  C- ���A�h���  C,�ʿ�p��h�v��  C,�:�����h���   C-�ْ:�h�x�@  C,~ߣ5��hi��  C,L��Z��hRz��  C,E���&�h:��   C,b����'�h#|�@  C,�md��h�ր  C,��T4�g�~��  C-R8�p��g���   C-�.u�gŀ�@  C-m"��pZ�g�ǀ  C-m8�0��g����  C,�� Fʽ�g�   C,��f: �gg��@  C,��7��gP��  C,�P��?�g8���  C,���z�g!�   C,�^�W^�g	��@  C,hI+i��f�	��  C,W��%���fڊ��  C,m�+�;R�f��   C,�/�;���f���@  C-������f���  C-�Q$ۃ�f|���  C.��)�)�fe�   C-�T5�5S�fM��@  C-������f6��  C-���2���f���  C-��:\�.�f�   C-��D��e@  C.
Bt�2W�e�|�  C.JY��Ϫ�e��x�  C-�IΥ>u�e�u   C-����Np�e��q@  C-��i,�ezm�  C.%J�1k�eb�i�  C.�tjB��eKf   C/�V'���e3�b@  C/�����e^�  C0�)[��e�Z�  C0��Q��d�W   C/�7����dՠS@  C0	�9`I=�d�!O�  C0`���d��K�  C0N�Z֊�d�#H   C/��xv���dw�D@  C.��3��d`%@�  C.U���dH�<�  C-�V]����d1'9   C-��U>!��d�5@  C-��H3�d)1�  C-�3j��)�c�-�  C-�d����c�+*   C-ج"UW�c��&@  C-��3G%��c�-"�  C-�\f�D��c���  C-?���2!�cu/   C-��W>�c]�@  C,���=x��cF1�  C,�2�u��c.��  C,��4�F��c3   C,\��6F��b��@  C,

�]��b�5�  C+���^$M�bж �  C,5�i���b�6�   C,:�ތ��b���@  C,�P��b�8��  C-8̏���br���  C-�-��b[:�   C-����m�bC��@  C-��[SfC�b,<�  C-�-�c5�b���  C-�7)��b�a�>�   C-� Z5�a��@  C.'.�	�b�a�@׀  C.g�&^"�a����  C.�m����a�B�   C.a�d����a���@  C.Hż��7�apDȀ  C.\��	�|�aX���  C.����;�aAF�   C/�U`=T�a)ǽ@  C/��ek���aH��  C0 �@�P��`�ɵ�  C0b+��H��`�J�   C0i��=��`�ˮ@  C0\��"�d�`�L��  C0m��BS�`�ͦ�  C0|f�!j%�`�N�   C0�*�����`mϟ@  C0v��"W8�`VP��  C0\�(2:�`>ї�  C//*U,s��`'R�   C.�-�P7M�`Ӑ@  C."4'Z��_�   C-��\��_���  C-�=w?|��_��
   C-�Eyᲊ�_c��  C-��v0�_4��   C.�z-B�_��  C-������^ִ�   C-<e����^���  C-��q��^x��   C,�Ż�w&�^I�Հ  C,�=B��<�^��   C,�"���]�ƀ  C,'JÊ��]���   C,�����]�·�  C+�g����]^İ   C,ċ0�]�]/ƨ�  C,3N�����] ȡ   C,y��M��\�ʙ�  C,��2�)��\�̒   C-��A�Ԓ�\sΊ�  C.�F�\DЃ   C.yQ�a�\�{�  C-��{����[��t   C-��p��[��l�  C-�N@U��[��e   C.��CA�[Y�]�  C.J�L���[*�V   C.b)ڑa$�Z��N�  C.�<o�Z��G   C-������Z��?�  C.��6N�Zn�8   C.T�Y�Z{�Z?�0�  C.�n(��[�Z�)   C/8`3��Y��!�  C/��*I/��Y��   C0%v)؈!�Y���  C0!}z�-�YT�   C0	�gy��Y%��  C0�n�_�X���   C0)ve����X���  C06�9rV\�X���   C0[֦	���Xi��  C/��Y����X:��   C.�����,�X�ր  C.PO�kAM�W���   C-�����T�W�ǀ  C-�$Rvr�W�   C-�ίf�3�WP��  C-�VA�Vo�W!�   C-��K`��V�	��  C-錙��V��   C-��4��6�V���  C-!��~�.�Ve�   C,�k�#^�V6��  C,�u�u�.�V�   C,��b�b�U�|�  C,��§u��U�u   C,�i_���Uzm�  C,,���UKf   C,��<��U^�  C,��aU	�T�W   C,Mu}�!�T�!O�  C,��H�3�T�#H   C,�
�(��T`%@�  C-T෡���T1'9   C-����d��T)1�  C-�9���5�S�+*   C-�>(7���S�-"�  C-��{�S'�Su/   C-�L6xU�SF1�  C.$M�����S3   C.k��@<��R�5�  C.� ��U��R�6�   C.�L��_�R�8��  C.���Y���R[:�   C.��*�|�R,<�  C.�$Yһ��Q�>�   C/YN���Q�@׀  C/�7���2�Q�B�   C0:��)�E�QpDȀ  C0vo侩~�QAF�   C0}�����QH��  C0oj��,��P�J�   C0}�T��y�P�L��  C0�gR,�P�N�   C0��v>n��PVP��  C0�����P'R�   C0\-8�6n�O�   C0WY�q>B�O��
   C0`��{��O4��   C0p:`���Nִ�   C0� �b���Nx��   C0�;�es�N��   C0쫴�f��M���   C112;�^B�M^İ   C1e��L��M ȡ   C1`;��b�L�̒   C1B�>�:��LDЃ   C1?�OS�K��t   C0�F	�2��K��e   C/�{_��z�K*�V   C-�3���J��G   C-�=�D�Jn�8   C,�5�!���J�)   C,9�Tb3�I��   C,&|�L�w�IT�   C,7�v��H���   C,m*n����H���   C,�m��~*�H:��   C-��e�A�G���   C-�U����G�   C-�F�1f��G!�   C-Ri�U�F��   C-Q����Fe�   C-`�M-�Q�F�   C-�N�l�E�u   C-�z�nN�EKf   C-���K��D�W   C-�������D�#H   C-�G�s2�D1'9   C-��al���C�+*   C-�Q����Cu/   C.>�N#�.�C3   C.��r��B�6�   C/]�����B[:�   C/ۨ�r��A�>�   C/Ԥ����A�B�   C/��+����AAF�   C/�@'N�@�J�   C/Ո�0O�@�N�   C0�1�@'R�   C0Lm"�x��?��
   C03�YJm �>ִ�   C02�>��>��   C0A,jB}�=^İ   C0XS��N�<�̒   C0y��1��;��t   C0��F���;*�V   C0������:n�8   C1.[�����9��   C1e�D����8���   C1j�4�2�8:��   C1V��&o��7�   C1UHOp���6��   C1W�;�.��6�   C1c�qM��5Kf   C1&���4�#H   C1��s��3�+*   C1���z�33   C1��k���2[:�   C1+W	��1�B�   C1LZ�P��0�J�   C1y�ʭ�w�0'R�   C1������.ִ�   C1��u;�-^İ   C23�Y���+��t   C25�L,�:�*n�8   C2!_�g�(���   C1`X&_Uy�'�   C/\ó���&�   C-����7�$�#H   C,�<
���#3   C,N�u�r�!�B�   C+� `"�� 'R�   C+�)�����^İ   C+�0|��n�8   C,;+bX��   C,�j�� ��#H   C-`i͙���B�   C-��4�W�^İ   C.6HdTL��   C.C��2!���B�   C.�U�����   C-��sX�A���   C. ��UC�        C.��.A��   C.J��]$�A��   C-����B�B�   C-���fBB�   C-���S
NB^İ   C. �"��B�B�   C.4�x���B�#H   C.� s�[B�   C/
&Yе�Bn�8   C/��R"oB^İ   C0
m�@�B 'R�   C0&ʪ�B!�B�   C/ޖR��	B#3   C/�;9s�zB$�#H   C0�l��>B&�   C0)��2�MB'�   C0f�#�YB(���   C/�Jh�B*n�8   C.�Bx/��B+��t   C.dc�� PB-^İ   C-�V�JRB.ִ�   C-�c�C�B0'R�   C-����0B0�J�   C-�o�ql�B1�B�   C.4�6��B2[:�   C.Vd���0B33   C.�]��B3�+*   C-�nA�B4�#H   C-o{��5oB5Kf   C-N���e�B6�   C-1��� B6��   C-)�>��B7�   C-#��G�zB8:��   C,���I��B8���   C,XP��"B9��   C,W,<S�oB:n�8   C,o~�bOB;*�V   C,���rPB;��t   C,����PMB<�̒   C-^il���B=^İ   C-�MP��B>��   C-��|�B>ִ�   C-=1NB?��
   C-%�_��B@'R�   C-tw���B@�N�   C-`����B@�J�   C--N�=BAAF�   C-G/���BA�B�   C,�s�w|�BA�>�   C,����9BB[:�   C,�����BB�6�   C,԰���EBC3   C-+:�a�DBCu/   C-�p7�}BC�+*   C.��`��BD1'9   C.s� <��BD�#H   C.k,;HƏBD�W   C.4�!�E�BEKf   C.8�u＊BE�u   C.D�{�YHBF�   C.b��Z�BFe�   C.�p�K߆BF��   C.�>=��fBG!�   C._`�_�wBG�   C.M�Ƀ�BG���   C.G�f���BH:��   C.L�"</�BH���   C.|�Em�BH���   C.��ahBIT�   C/ �X^`BI��   C/)�BՍNBJ�)   C.ԙ��Z�BJn�8   C.>��B�BJ��G   C-��T���BK*�V   C-�`3���BK��e   C-r6�$�BK��t   C-^rUи�BLDЃ   C-!4���BL�̒   C,��X�1BM ȡ   C,Y�:P��BM^İ   C,M�j�ߡBM���   C,i � �BN��   C,�U��xBNx��   C,����J�BNִ�   C-6+��BO4��   C-s�yN�BO��
   C-GS6��BO�   C,�)b	�hBP'R�   C,�<A�Z�BPVP��  C,�VLe:XBP�N�   C,��ևSBP�L��  C,�[IBP�J�   C,�~we�BQH��  C,x���0BQAF�   C,1�2�JBQpDȀ  C,,T��+ BQ�B�   C,T�G�uBQ�@׀  C,��$��RBQ�>�   C-��*N6BR,<�  C-������BR[:�   C-�/��zBR�8��  C-�\=ݱBR�6�   C-�x}'��BR�5�  C-�-e�!�BS3   C-������BSF1�  C.#l�BSu/   C.[.�_ձBS�-"�  C.�o�#�IBS�+*   C.Z�v�ڳBT)1�  C.A=�24�BT1'9   C.ZK��BT`%@�  C.����BT�#H   C/�w[q�BT�!O�  C/�g��'�BT�W   C0��i�BU^�  C0Q��8�BUKf   C0Uz�R1BUzm�  C03(��BU�u   C0*�1��BU�|�  C0*V q�cBV�   C0/�WŽ�BV6��  C0 �nG�rBVe�   C/�z�� �BV���  C/گ��fqBV��   C/��b{wBV�	��  C0	&���GBW!�   C0+���ǆBWP��  C0Z����BW�   C0���e0BW�ǀ  C0�,�_�BW���   C1�-�R�BX�ր  C1�&��PBX:��   C0׺\@BXi��  C/�Ӹ��BX���   C/<3���BX���  C-�s�v�3BX���   C,ǝ�0�BY%��  C,&R��wSBYT�   C+�I�s�BY���  C+��j�XBY��   C+})	Y1BY��!�  C+�(:4�BZ�)   C+�2��R�BZ?�0�  C,_׶�:�BZn�8   C-
�U���BZ��?�  C-�ٟ�̸BZ��G   C-��R�h�BZ��N�  C-V=
r�B[*�V   C-X0���KB[Y�]�  C-[�T(B[��e   C-��' 2�B[��l�  C-����MB[��t   C-����ΎB\�{�  C.
�A�n1B\DЃ   C-͖��C B\sΊ�  C-�kzY�B\�̒   C-�1����B\�ʙ�  C.\�f~�HB] ȡ   C.�K���B]/ƨ�  C/{0봈�B]^İ   C/�"��B]�·�  C/�m!�|B]���   C/�����.B]�ƀ  C/��B-@B^��   C/���dB^I�Հ  C/ݨ;�$B^x��   C0$��O��B^���  C0�O�#B^ִ�   C0
�?�}�B_��  C0�z�)�B_4��   C03@0Q�B_c��  C0Y�KEB_��
   C0�痰h*B_���  C0�KΣ�B_�   C1�����B`Ӑ@  C1ESx�hB`'R�   C1R%��WUB`>ї�  C1Ip:�h?B`VP��  C1YSߴv%B`mϟ@  C1fP|O�IB`�N�   C1u�	\u�B`�ͦ�  C1p��%�B`�L��  C1:9�FhB`�ˮ@  C1*�c�/B`�J�   C0>�����B`�ɵ�  C/0x����BaH��  C.f�����Ba)ǽ@  C-�
���oBaAF�   C-�`�D׹BaX���  C-鶃rr�BapDȀ  C-�t�)RBa���@  C-�y����Ba�B�   C-5-;֬Ba����  C-	-UB�Ba�@׀  C,�,/�Ba��@  C,؛�m��Ba�>�   C,��[��{Bb���  C,��"�Bb,<�  C,I�7d�;BbC��@  C,%�Y,Bb[:�   C,,�=}��Bbr���  C,\>Xw�Bb�8��  C,��<�d�Bb���@  C-5`����Bb�6�   C-�7�Z��Bbж �  C.S��h�Bb�5�  C.O}z��=Bb��@  C. ��wIBc3   C-�-�6WBc.��  C-݆-�hBcF1�  C-��k�SBc]�@  C.&�6��Bcu/   C.`�W�*7Bc���  C.Q,� UiBc�-"�  C.ow��Bc��&@  C.R&��tBc�+*   C.>�b3yBc�-�  C.���u��Bd)1�  C/ ����Bd�5@  C/�oFl�Bd1'9   C0ɕ�>HBdH�<�  C0��	Bd`%@�  C0 �H�LBdw�D@  C0�-8�Bd�#H   C0��]�JBd��K�  C0.��{E{Bd�!O�  C/��µf_BdՠS@  C.͡�=�Bd�W   C.�(D1�Be�Z�  C-�1U��Be^�  C->�B�Be3�b@  C- ���BeKf   C-%��GV�Beb�i�  C-P$��fBezm�  C-�"�2<�Be��q@  C-ԣu<�Be�u   C-�9H'�Be��x�  C-Q=��}Be�|�  C--V4e��Be@  C-�FRBf�   C-ʝ'�Bf���  C-i\� 6Bf6��  C-.� ���BfM��@  C,�h��tBfe�   C,v��#�Bf|���  C,up�d�Bf���  C,�w&)�Bf���@  C,��HGN{Bf��   C-or �Z�Bfڊ��  C-��W�P�Bf�	��  C.V�b>��Bg	��@  C.Q ,��xBg!�   C.
��e[Bg8���  C.�Z���BgP��  C.
.`_�Bgg��@  C.>�0��Bg�   C.|�.w�iBg����  C.����� Bg�ǀ  C.t�ʇ�KBgŀ�@  C.fh!Bg���   C.�X�?E0Bg�~��  C.�V���]Bh�ր  C/?6� ��Bh#|�@  C/̈́���Bh:��   C01���BhRz��  C0o���}Bhi��  C0x\���Bh�x�@  C0fV����Bh���   C0jz�K�Bh�v��  C0q�R�Bh���  C0��m��Bh�t�@  C0q�ɢLBh���   C0F�&��Bir��  C0?D����Bi%��  C0DC2��CBi=q@  C0R>z�ӶBiT�   C0p�H<��Bilo�  C0��T*�mBi���  C0܇�eF�Bi�m@  C1#��hcBi��   C1]f���Bi�k�  C1[��*�nBi��!�  C1?�n��Bi�i%@  C0H�=�T�Bj�)   C/y�o�i�Bj(g,�  C.�i��c�Bj?�0�  C-;�`S8�BjWe4@  C,��q��0Bjn�8   C,,B��Bj�c;�  C,3g�BQBj��?�  C,��ث�Bj�aC@  C,K}���[Bj��G   C,��@+�#Bj�_J�  C,����Bj��N�  C-~�I�Bk]R@  C. dw(g�Bk*�V   C-�J��BkB[Y�  C-�\|�I�BkY�]�  C-����BkqYa@  C-~�%��Bk��e   C-���G��Bk�Wh�  C-����JBk��l�  C-�r����Bk�Up@  C-�����qBk��t   C-��z�̠Bk�Sw�  C-��{�Bl�{�  C-��W�BMBl-Q@  C.LU۝�BlDЃ   C.�Sj��jBl\O��  C/j�W`��BlsΊ�  C/�>뢳.Bl�M�@  C/����:Bl�̒   C/��u���Bl�K��  C/�7"� hBl�ʙ�  C/���eӯBl�I�@  C0�!Bm ȡ   C000�r*�BmG��  C/6b�:Bm/ƨ�  C.�a�OBmGE�@  C.	94���Bm^İ   C-�V��BmvC��  C-��(�5�Bm�·�  C-���eGBm�A�@  C-�%���mBm���   C.'�d^�Bm�?��  C.K�+�/Bm�ƀ  C.
�K�FMBn=�@  C-�)pf��Bn��   C-g:
J9_Bn2;��  C-;@��\�BnI�Հ  C-(��:\Bna9�@  C-|̟�,Bnx��   C,�3�Ӱ�Bn�7��  C,[V)�Bn���  C,?.@N�FBn�5�@  C,:'qX�Bnִ�   C,ZKo��Bn�3��  C,�p�\OBo��  C-�����Bo1�@  C-w��3v�Bo4��   C-�dN2#BoL/��  C-����!PBoc��  C-��%\l�Bo{.@  C-��OENBo��
   C-x�ܛBo�,�  C-����*Bo���  C-�A�4=vBo�*@  C.�)��Bo�   C-��o:[�Bp`  C-���[�BpӐ@  C-�A��fX