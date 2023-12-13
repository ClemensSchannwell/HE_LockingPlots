CDF  �   
      time       lon       lat          	   CDI       @Climate Data Interface version 2.0.3 (https://mpimet.mpg.de/cdi)   Conventions       CF-1.5     source        #PISM development no-version-control    command      / /home/m/m300792/ModelCodes/mPISM_build/install_03-10-22/bin/pismr -i pism_-065000/pism_-065000.nc -stress_balance ssa+sia -y 100 -ys -65000 -o pism_-064900/pism_-064900.nc -surface_given_file forcing/smb_-064900.nc -surface given -ocean_th_file forcing/oce_-064900.nc -ocean th -topg_file pism_-065000/pism_-065000_new_topg.nc -bed_def prescribed -uplift_file pism_-065000/pism_-065000_dbdt.nc -extra_times 10 -extra_vars thk,usurf,bheatflx,velsurf_mag,velbase_mag,tillwat,bmelt,velbase,flux_divergence,velsurf,velbar,wvelsurf,wvelbase,climatic_mass_balance_cumulative,potential_climatic_mass_balance_cumulative,tempbase,tempsurf,mask,temppabase,tempicethk_basal,topg,nonneg_flux_cumulative,grounded_basal_flux_cumulative,floating_basal_flux_cumulative,discharge_flux_cumulative,Href,taud_mag,lat,lon,calving_mask,sea_level -extra_file pism_-064900/pism_-064900_extra -extra_split -ts_times yearly -ts_file pism_-064900/pism_-064900_ts.nc -topg_to_phi 15.0,30.0,-300.0,400.0 -ssafd_ksp_rtol 1e-7 -enthalpy_temperate_conductivity_ratio 0.001 -pseudo_plastic -pseudo_plastic_q .25 -pseudo_plastic_uthreshold 70 -till_effective_fraction_overburden 0.04 -tauc_slippery_grounding_lines -sia_e 8 -ssa_e 1 -ssa_upwind_factor 0.7 -calving eigen_calving,thickness_calving -thickness_calving_threshold 100 -part_redist -part_grid -cfbc -kill_icebergs -eigen_calving_K 1e16 -verbose 2 -sliding_scale_factor_reduces_tauc 1 -th_heat_coefficient 4e-5 -th_salinity_coefficient 4e-8 -o_format netcdf3 -sediment_file till_mask.nc -no_divide_tillphi -periodicity none -ocean_th_period 1 -o_order yxz
    proj4         j+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs     history_of_appended_files        �Fri Jul 29 07:48:43 2022: Appended file outFinal2 had following "history" attribute:
Fri Jul 29 07:48:41 2022: ncap2 -O -s lon_bnds=float(lon_bnds) outFinal1 outFinal2
Fri Jul 29 07:48:41 2022: ncap2 -O -s lat_bnds=float(lat_bnds) outFinal outFinal1
Fri Jul 29 07:48:40 2022: cdo -add thkNew thkOld outFinal
Fri Jul 29 07:48:38 2022: cdo -mul -selvar,thk pism_-064900/pism_-064900.nc MaskLaurentide thkNew
     NCO       _netCDF Operators version 5.0.6 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)    history      7�Wed Sep 14 12:10:35 2022: cdo -mergetime File000001.nc File000002.nc File000003.nc File000004.nc File000005.nc File000006.nc File000007.nc File000008.nc File000009.nc File000010.nc File000011.nc File000012.nc File000013.nc File000014.nc File000015.nc File000016.nc File000017.nc File000018.nc File000019.nc File000020.nc File000021.nc File000022.nc File000023.nc File000024.nc File000025.nc File000026.nc File000027.nc File000028.nc File000029.nc File000030.nc File000031.nc File000032.nc File000033.nc File000034.nc File000035.nc File000036.nc File000037.nc File000038.nc File000039.nc File000040.nc File000041.nc File000042.nc File000043.nc File000044.nc File000045.nc File000046.nc File000047.nc File000048.nc File000049.nc File000050.nc File000051.nc File000052.nc File000053.nc File000054.nc File000055.nc File000056.nc File000057.nc File000058.nc File000059.nc File000060.nc File000061.nc File000062.nc File000063.nc File000064.nc File000065.nc File000066.nc File000067.nc File000068.nc File000069.nc File000070.nc File000071.nc File000072.nc File000073.nc File000074.nc File000075.nc File000076.nc File000077.nc File000078.nc File000079.nc File000080.nc File000081.nc File000082.nc File000083.nc File000084.nc File000085.nc File000086.nc File000087.nc File000088.nc File000089.nc File000090.nc File000091.nc File000092.nc File000093.nc File000094.nc File000095.nc File000096.nc File000097.nc File000098.nc File000099.nc File000100.nc File000101.nc File000102.nc File000103.nc File000104.nc File000105.nc File000106.nc File000107.nc File000108.nc File000109.nc File000110.nc File000111.nc File000112.nc File000113.nc File000114.nc File000115.nc File000116.nc File000117.nc File000118.nc File000119.nc File000120.nc File000121.nc File000122.nc File000123.nc File000124.nc File000125.nc File000126.nc File000127.nc File000128.nc File000129.nc File000130.nc File000131.nc File000132.nc File000133.nc File000134.nc File000135.nc File000136.nc File000137.nc File000138.nc File000139.nc File000140.nc File000141.nc File000142.nc File000143.nc File000144.nc File000145.nc File000146.nc File000147.nc File000148.nc File000149.nc File000150.nc File000151.nc File000152.nc File000153.nc File000154.nc File000155.nc File000156.nc File000157.nc File000158.nc File000159.nc File000160.nc File000161.nc File000162.nc File000163.nc File000164.nc File000165.nc File000166.nc File000167.nc File000168.nc File000169.nc File000170.nc File000171.nc File000172.nc File000173.nc File000174.nc File000175.nc File000176.nc File000177.nc File000178.nc File000179.nc File000180.nc File000181.nc File000182.nc File000183.nc File000184.nc File000185.nc File000186.nc File000187.nc File000188.nc File000189.nc File000190.nc File000191.nc File000192.nc File000193.nc File000194.nc File000195.nc File000196.nc File000197.nc File000198.nc File000199.nc File000200.nc File000201.nc File000202.nc File000203.nc File000204.nc File000205.nc File000206.nc File000207.nc File000208.nc File000209.nc File000210.nc File000211.nc File000212.nc File000213.nc File000214.nc File000215.nc File000216.nc File000217.nc File000218.nc File000219.nc File000220.nc File000221.nc File000222.nc File000223.nc File000224.nc File000225.nc File000226.nc File000227.nc File000228.nc File000229.nc File000230.nc File000231.nc File000232.nc File000233.nc File000234.nc File000235.nc File000236.nc File000237.nc File000238.nc File000239.nc File000240.nc File000241.nc File000242.nc File000243.nc File000244.nc File000245.nc File000246.nc File000247.nc File000248.nc File000249.nc File000250.nc File000251.nc File000252.nc File000253.nc File000254.nc File000255.nc File000256.nc File000257.nc File000258.nc File000259.nc File000260.nc File000261.nc File000262.nc File000263.nc File000264.nc File000265.nc File000266.nc File000267.nc File000268.nc File000269.nc File000270.nc File000271.nc File000272.nc File000273.nc File000274.nc File000275.nc File000276.nc File000277.nc File000278.nc File000279.nc File000280.nc File000281.nc File000282.nc File000283.nc File000284.nc File000285.nc File000286.nc File000287.nc File000288.nc File000289.nc File000290.nc File000291.nc File000292.nc File000293.nc File000294.nc File000295.nc File000296.nc File000297.nc File000298.nc File000299.nc File000300.nc File000301.nc File000302.nc File000303.nc File000304.nc File000305.nc File000306.nc File000307.nc File000308.nc File000309.nc File000310.nc File000311.nc File000312.nc File000313.nc File000314.nc File000315.nc File000316.nc File000317.nc File000318.nc File000319.nc File000320.nc File000321.nc File000322.nc File000323.nc File000324.nc File000325.nc File000326.nc File000327.nc File000328.nc File000329.nc File000330.nc File000331.nc File000332.nc File000333.nc File000334.nc File000335.nc File000336.nc File000337.nc File000338.nc File000339.nc File000340.nc File000341.nc File000342.nc File000343.nc File000344.nc File000345.nc File000346.nc File000347.nc File000348.nc File000349.nc File000350.nc File000351.nc File000352.nc File000353.nc File000354.nc File000355.nc File000356.nc File000357.nc File000358.nc File000359.nc File000360.nc File000361.nc File000362.nc File000363.nc File000364.nc File000365.nc File000366.nc File000367.nc File000368.nc File000369.nc File000370.nc File000371.nc File000372.nc File000373.nc File000374.nc File000375.nc File000376.nc File000377.nc File000378.nc File000379.nc File000380.nc File000381.nc File000382.nc File000383.nc File000384.nc File000385.nc File000386.nc File000387.nc File000388.nc File000389.nc File000390.nc File000391.nc File000392.nc File000393.nc File000394.nc File000395.nc File000396.nc File000397.nc File000398.nc File000399.nc File000400.nc File000401.nc File000402.nc File000403.nc File000404.nc File000405.nc File000406.nc File000407.nc File000408.nc File000409.nc File000410.nc File000411.nc File000412.nc File000413.nc File000414.nc File000415.nc File000416.nc File000417.nc File000418.nc File000419.nc File000420.nc File000421.nc File000422.nc File000423.nc File000424.nc File000425.nc File000426.nc File000427.nc File000428.nc File000429.nc File000430.nc File000431.nc File000432.nc File000433.nc File000434.nc File000435.nc File000436.nc File000437.nc File000438.nc File000439.nc File000440.nc File000441.nc File000442.nc File000443.nc File000444.nc File000445.nc File000446.nc File000447.nc File000448.nc File000449.nc File000450.nc File000451.nc File000452.nc File000453.nc File000454.nc File000455.nc File000456.nc File000457.nc File000458.nc File000459.nc File000460.nc File000461.nc File000462.nc File000463.nc File000464.nc File000465.nc File000466.nc File000467.nc File000468.nc File000469.nc File000470.nc File000471.nc File000472.nc File000473.nc File000474.nc File000475.nc File000476.nc File000477.nc File000478.nc File000479.nc File000480.nc File000481.nc File000482.nc File000483.nc File000484.nc File000485.nc File000486.nc File000487.nc File000488.nc File000489.nc File000490.nc File000491.nc File000492.nc File000493.nc File000494.nc File000495.nc File000496.nc File000497.nc File000498.nc File000499.nc File000500.nc File000501.nc File000502.nc File000503.nc File000504.nc File000505.nc File000506.nc File000507.nc File000508.nc File000509.nc File000510.nc File000511.nc File000512.nc File000513.nc File000514.nc File000515.nc File000516.nc File000517.nc File000518.nc File000519.nc File000520.nc File000521.nc File000522.nc File000523.nc File000524.nc File000525.nc File000526.nc File000527.nc File000528.nc File000529.nc File000530.nc File000531.nc File000532.nc File000533.nc File000534.nc File000535.nc File000536.nc File000537.nc File000538.nc File000539.nc File000540.nc File000541.nc File000542.nc File000543.nc File000544.nc File000545.nc File000546.nc File000547.nc File000548.nc File000549.nc File000550.nc File000551.nc File000552.nc File000553.nc File000554.nc File000555.nc File000556.nc File000557.nc File000558.nc File000559.nc File000560.nc File000561.nc File000562.nc File000563.nc File000564.nc File000565.nc File000566.nc File000567.nc File000568.nc File000569.nc File000570.nc File000571.nc File000572.nc File000573.nc File000574.nc File000575.nc File000576.nc File000577.nc File000578.nc File000579.nc File000580.nc File000581.nc File000582.nc File000583.nc File000584.nc File000585.nc File000586.nc File000587.nc File000588.nc File000589.nc File000590.nc File000591.nc File000592.nc File000593.nc File000594.nc File000595.nc File000596.nc File000597.nc File000598.nc File000599.nc File000600.nc File000601.nc File000602.nc File000603.nc File000604.nc File000605.nc File000606.nc File000607.nc File000608.nc File000609.nc File000610.nc File000611.nc File000612.nc File000613.nc File000614.nc File000615.nc File000616.nc File000617.nc File000618.nc File000619.nc File000620.nc File000621.nc File000622.nc File000623.nc File000624.nc File000625.nc File000626.nc File000627.nc File000628.nc File000629.nc File000630.nc File000631.nc File000632.nc File000633.nc File000634.nc File000635.nc File000636.nc File000637.nc File000638.nc File000639.nc File000640.nc File000641.nc File000642.nc File000643.nc File000644.nc File000645.nc File000646.nc File000647.nc File000648.nc File000649.nc File000650.nc File000651.nc File000652.nc File000653.nc File000654.nc File000655.nc File000656.nc File000657.nc File000658.nc File000659.nc File000660.nc File000661.nc File000662.nc File000663.nc File000664.nc File000665.nc File000666.nc File000667.nc File000668.nc File000669.nc File000670.nc File000671.nc File000672.nc File000673.nc File000674.nc File000675.nc File000676.nc File000677.nc File000678.nc File000679.nc File000680.nc File000681.nc File000682.nc File000683.nc File000684.nc File000685.nc File000686.nc File000687.nc File000688.nc File000689.nc File000690.nc File000691.nc File000692.nc File000693.nc File000694.nc File000695.nc File000696.nc File000697.nc File000698.nc File000699.nc File000700.nc File000701.nc File000702.nc File000703.nc File000704.nc File000705.nc File000706.nc File000707.nc File000708.nc File000709.nc File000710.nc File000711.nc File000712.nc File000713.nc File000714.nc File000715.nc File000716.nc File000717.nc File000718.nc File000719.nc File000720.nc File000721.nc File000722.nc File000723.nc File000724.nc File000725.nc File000726.nc File000727.nc File000728.nc File000729.nc File000730.nc File000731.nc File000732.nc File000733.nc File000734.nc File000735.nc File000736.nc File000737.nc File000738.nc File000739.nc File000740.nc File000741.nc File000742.nc File000743.nc File000744.nc File000745.nc File000746.nc File000747.nc File000748.nc File000749.nc File000750.nc File000751.nc File000752.nc File000753.nc File000754.nc File000755.nc File000756.nc File000757.nc File000758.nc File000759.nc File000760.nc File000761.nc File000762.nc File000763.nc File000764.nc File000765.nc File000766.nc File000767.nc File000768.nc File000769.nc File000770.nc File000771.nc File000772.nc File000773.nc File000774.nc File000775.nc File000776.nc File000777.nc File000778.nc File000779.nc File000780.nc File000781.nc File000782.nc File000783.nc File000784.nc File000785.nc File000786.nc File000787.nc File000788.nc File000789.nc File000790.nc File000791.nc File000792.nc File000793.nc File000794.nc File000795.nc File000796.nc File000797.nc File000798.nc File000799.nc File000800.nc File000801.nc File000802.nc File000803.nc File000804.nc File000805.nc File000806.nc File000807.nc File000808.nc File000809.nc File000810.nc File000811.nc File000812.nc File000813.nc File000814.nc File000815.nc File000816.nc File000817.nc File000818.nc File000819.nc File000820.nc File000821.nc File000822.nc File000823.nc File000824.nc File000825.nc File000826.nc File000827.nc File000828.nc File000829.nc File000830.nc File000831.nc File000832.nc File000833.nc File000834.nc File000835.nc File000836.nc File000837.nc File000838.nc File000839.nc File000840.nc File000841.nc File000842.nc File000843.nc File000844.nc File000845.nc File000846.nc File000847.nc File000848.nc File000849.nc File000850.nc File000851.nc File000852.nc File000853.nc File000854.nc File000855.nc File000856.nc File000857.nc File000858.nc File000859.nc File000860.nc File000861.nc File000862.nc File000863.nc File000864.nc File000865.nc File000866.nc File000867.nc File000868.nc File000869.nc File000870.nc File000871.nc File000872.nc File000873.nc File000874.nc File000875.nc File000876.nc File000877.nc File000878.nc File000879.nc File000880.nc File000881.nc File000882.nc File000883.nc File000884.nc File000885.nc File000886.nc File000887.nc File000888.nc File000889.nc File000890.nc File000891.nc File000892.nc File000893.nc File000894.nc File000895.nc File000896.nc File000897.nc File000898.nc File000899.nc File000900.nc File000901.nc File000902.nc File000903.nc File000904.nc File000905.nc File000906.nc File000907.nc File000908.nc File000909.nc File000910.nc File000911.nc File000912.nc File000913.nc File000914.nc File000915.nc File000916.nc File000917.nc File000918.nc File000919.nc File000920.nc File000921.nc File000922.nc File000923.nc File000924.nc File000925.nc File000926.nc File000927.nc File000928.nc File000929.nc File000930.nc File000931.nc File000932.nc File000933.nc File000934.nc File000935.nc File000936.nc File000937.nc File000938.nc File000939.nc File000940.nc File000941.nc File000942.nc File000943.nc File000944.nc File000945.nc File000946.nc File000947.nc File000948.nc File000949.nc File000950.nc File000951.nc File000952.nc File000953.nc File000954.nc File000955.nc File000956.nc File000957.nc File000958.nc File000959.nc File000960.nc File000961.nc File000962.nc File000963.nc File000964.nc File000965.nc File000966.nc File000967.nc File000968.nc File000969.nc File000970.nc File000971.nc File000972.nc File000973.nc File000974.nc File000975.nc File000976.nc File000977.nc File000978.nc File000979.nc File000980.nc File000981.nc File000982.nc File000983.nc File000984.nc File000985.nc File000986.nc File000987.nc File000988.nc File000989.nc File000990.nc File000991.nc File000992.nc File000993.nc File000994.nc File000995.nc File000996.nc File000997.nc File000998.nc File000999.nc File001000.nc IceVolumeKenzie.nc
Wed Sep 14 12:04:15 2022: cdo -fldsum -mulc,10000 -mulc,10000 -selindexbox,211,427,600,415 -selvar,thk /work/ba0989/m300792/HE_Runs/ExperimentsComposite//HE91//pism_-064900/pism_-064900.nc Tmp/File000001.nc     CDO       @Climate Data Operators version 2.0.3 (https://mpimet.mpg.de/cdo)         time                standard_name         time   	long_name         time   units         seconds since 1-1-1    calendar      365_day    axis      T               D�   lon                standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X               D�   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y               D�   thk                       standard_name         land_ice_thickness     	long_name         land ice thickness     units         m      pism_intent       model_state             D�                �}Ȁ@�  C-V��z��}���   C-M�{���}�=   C->�9�V�}�A�@  C-4ڇU��}��9`  C-@��q�}�·�  C-a�]�e�}�5�  C-�#��6��}vC��  C-�t�����}j�1�  C-蒥��o�}^İ   C.+���u��}S.   C.r�>���}GE�@  C.�R�&��};�*`  C/��dA�}/ƨ�  C/�r���`�}$&�  C/�
A�Z�}G��  C/��ʝ��}�"�  C0�:}\��} ȡ   C0)R;��|�	   C03�4���|�I�@  C08����|݊`  C03�5U��|�ʙ�  C0T�7��4�|��  C0p�;�fT�|�K��  C0�_��r�|���  C0��z:���|�̒   C0�o����|�   C0�>�9��|�M�@  C1�	?��|�`  C1<�>�P�|sΊ�  C1gaE-�|h�  C1r�c����|\O��  C1������|P��  C1��^t��|DЃ   C1�@�
���|9   C1f0hm�U�|-Q@  C1Xw��(�|!��`  C1T�k�i��|�{�  C1O������|
��  C0к)��{�Sw�  C0(�����{���  C/e�N����{��t   C.xs@��>�{��   C.>.�1N�{�Up@  C-ިey�l�{Õ�`  C-�0�Y�{��l�  C-g���{��  C-.^W�g�{�Wh�  C,��'��{����  C,ZW����{��e   C,--I'�D�{}�   C,�6K�w�{qYa@  C,�2 ;�{e��`  C,O{̯�{Y�]�  C,-D�L���{N۠  C,c��Or��{B[Y�  C,���i9�{6���  C-�
��(�{*�V   C-w��Z3�{�   C-��F�z�{]R@  C.Fu�SP��{��`  C.\H��B��z��N�  C.rfMqC�z�̠  C.wџ�ww�z�_J�  C.��D�o��z؟��  C.�̆1"�z��G   C.�����z� �   C.��FO�l�z�aC@  C.���mG�z���`  C.��p����z��?�  C.���6�~�z�"��  C.� �:��z�c;�  C/�G\�zz���  C/j�УM�zn�8   C/�Y5N��zc$�   C0���*��zWe4@  C0W��e�zK��`  C0/(�����z?�0�  C0?��`���z4&��  C0T��LF��z(g,�  C0zzk�.V�z���  C0�jeǩ��z�)   C0�A�ǔ�z(�   C0��K��y�i%@  C0�Q�p�y���`  C0��}�P�y��!�  C0��kC��y�*��  C0����A�y�k�  C1'q�P��y����  C1Y2�`_+�y��   C1��>\^��y�,�   C1�`�=��y�m@  C1��̙]�y���`  C1�Lޚ��y���  C0N�H-1��yx.��  C/Q'�a[�ylo�  C-�� �2��y`���  C,��CȋJ�yT�   C,� X��yI0�   C+�ݘ��?�y=q@  C+ۮ(�Al�y1��`  C+�����y%��  C,9��Q9��y2��  C,�Bc�x�yr��  C-���w�y�}�  C-v�fĶ!�x���   C-���<#��x�4z   C-�}�lU�x�t�@  C-��Y���xӵv`  C-�^z����x���  C-�t�L�x�6r�  C.���9�x�v��  C.I�sxQC�x��n�  C.i��y�x���   C.6��》�x�8k   C.D��*��x�x�@  C.b�U,��xu�g`  C.��c��4�xi��  C.��d��>�x^:c�  C/[��R��xRz��  C/�)`5���xF�_�  C/ȎQ-�x:��   C/��u�w<�x/<\   C0ڸS�R�x#|�@  C0�r9"-�x�X`  C0K�/k~/�x�ր  C0O[؈�x >T�  C0OP^�'�w�~��  C0XG٤��w�P�  C0fH�@3�w���   C0zڏB�w�@M   C0�,�v��wŀ�@  C0�U�G�N�w��I`  C0�X�e O�w�ǀ  C1�_�֏�w�BE�  C18��I�@�w����  C1FH�����w��A�  C1X�����w�   C1l��wsD>   C1yN��[�wg��@  C1���4��w[�:`  C1[6����wP��  C1Tqk���wDF6�  C1V��;���w8���  C1_���W��w,�2�  C1n5�w!�   C1�m�c���wH/   C1�`w�Z�w	��@  C1͂���B�v��+`  C1��-� �v�	��  C2$�d�;��v�J'�  C2,��$�I�vڊ��  C2;��_�E�v��#�  C1y�#���v��   C/�;�G4�v�L    C.�[��v���@  C-+�i�v��`  C,|�ݛ$��v���  C,4Ox�y*�v�N�  C,L����v|���  C,�>)R��vp��  C,�Ov�>�ve�   C-*'�X/�vYP   C-��.J���vM��@  C-���5@y�vA�`  C._�s?	��v6��  C.v+e�{�v*R	�  C.��.4�;�v���  C.��F�,(�v��  C.��w*�Y�v�   C.���a��u�T   C/	_,_}��u@  C/,c�V9��u���`  C.��G^A��u�|�  C.��C��u�U��  C.�mkMe��u��x�  C.�� ]��u����  C/;��z��u�u   C/�f�����u�W�   C/��+!��u��q@  C0(T�z� �u���`  C02�McP��uzm�  C0G��ʌ�unY�  C0b�K�ub�i�  C0h�.�p�uV���  C0�w��(�uKf   C0��u?r��u?[�   C0�[�*�J�u3�b@  C0�*�Qa��u'��`  C0��n��u^�  C0���D�u]ܠ  C0��]���u�Z�  C0���]��t����  C1/�h=��t�W   C1>����i�t�_�   C1i��K;��tՠS@  C1v|P6��t���`  C1����P�t�!O�  C1��z�O��t�a͠  C1������t��K�  C1��5z;�t����  C1����6��t�#H   C1���[��t�c�   C1��	�ĕ�tw�D@  C1�� �.�tk��`  C1���ʒ��t`%@�  C1�h|�?�tTe��  C1�K�d?��tH�<�  C1�j��|�t<��  C2���G�t1'9   C28HڿJ)�t%g�   C2?%����t�5@  C1}"��2��t�`  C04,�7�3�t)1�  C/��Y�"�s�i��  C-ֵ�D��s�-�  C,�
�-A��s���  C,p	�g�s�+*   C,>�9b�s�k�   C,@|F�s��&@  C,o��'���s��`  C,��٬���s�-"�  C-8���s�m��  C-dYVơ��s���  C-�O<A�]�s���  C.8�Mɏ�su/   C.N��曋�sio�   C.YK��s]�@  C.o[q5�F�sQ�`  C.�Azq��sF1�  C.��Y�Lm�s:q��  C.�c!r&��s.��  C/b�y���s"��  C.������s3   C.����Sp�ss�   C.�tW�*Z�r��@  C.�C�0�r��`  C/e���d�r�5�  C/fG���}�r�u��  C/�W�1M�rж �  C0(���i�r��~�  C0'Nq��r�6�   C00E��Wy�r�w{   C0I���%�r���@  C0Z��:��r��w`  C0��ɾ��r�8��  C0�){E���r~ys�  C0�rsE+�rr���  C0����!��rf�o�  C0�n����r[:�   C0�m�߮5�rO{l   C0҇,��W�rC��@  C0�o���r7�h`  C16�SF�r,<�  C1;��Uu�r }d�  C1f�u����r���  C1t)Q�9f�r�`�  C1��u+��q�>�   C1���=�
�q�]   C1���f��q��@  C1Ô4��q� Y`  C1�X�A4�q�@׀  C1~[�|ԇ�qU�  C1~��o��q����  C1�SI�=��q�Q�  C0�u-̼��q�B�   C/�l���q��N   C.�=�8��q���@  C.v���q|J`  C-���Eh�qpDȀ  C-���_D�qd�F�  C-Z�h�k�qX���  C-�R;��qMB�  C,�����qAF�   C,��=�W��q5�?   C,����� �q)ǽ@  C-
ɮ�	R�q;`  C-
��#�qH��  C,�%i���q�7�  C,��J�x"�p�ɵ�  C,���IV�p�
3�  C,������p�J�   C-Az�W�p׋0   C-����#�p�ˮ@  C.J�����p�,`  C._��@���p�L��  C.t:��p��(�  C.���:��p�ͦ�  C.�1$���p�$�  C.��� `�p�N�   C/}VZ���py�!   C/EI7���pmϟ@  C/~r5kb�pb`  C/e	�����pVP��  C//x���pJ��  C/6����P�p>ї�  C/U���p3�  C/�<��#9�p'R�   C/�Gik��p�   C0+��_"�pӐ@  C0YRv��[�p`  C0g:�r2�o�   C0��k,�o�*@  C0������o���  C0�obf|��o�,�  C0ט�]���o��
   C0�¤�O�o{.@  C0�1��0�oc��  C0ǐL�V��oL/��  C0�gl(4)�o4��   C0� �8��o1�@  C0"tz"MD�o��  C/-�m���n�3��  C.�(����nִ�   C.E	��a��n�5�@  C-�E��n���  C-���b��n�7��  C-<�]`�.�nx��   C-^|a���na9�@  C,�Kh��b�nI�Հ  C,צ�����n2;��  C,�7bd�n�n��   C,�e�е�n=�@  C,b.���m�ƀ  C,9xC!���m�?��  C,KU�W�m���   C,�8�B���m�A�@  C,�M���m�·�  C-Ud׏RU�mvC��  C-��5)�m^İ   C.%��ٿj�mGE�@  C.8�A��m/ƨ�  C.Ey�a���mG��  C.^r�.B�m ȡ   C.{��?�<�l�I�@  C.�S�Jۢ�l�ʙ�  C.�������l�K��  C/TU�p��l�̒   C.�_>��l�M�@  C.�pf���lsΊ�  C.��Xj�l\O��  C.�ʞ��T�lDЃ   C.���@Ͱ�l-Q@  C/"���b�l�{�  C/�m�����k�Sw�  C/�n���k��t   C/��ߚ��k�Up@  C0k5�8�k��l�  C0�~��k�Wh�  C0=3�j��k��e   C0m&�{��kqYa@  C0���>��kY�]�  C0x"mT��kB[Y�  C0{^� L�k*�V   C0�=u@=��k]R@  C0� �#p�j��N�  C0��L���j�_J�  C0ك"����j��G   C1��.��j�aC@  C1;4�^��j��?�  C1h��:Id�j�c;�  C1{�JA��jn�8   C1R���O=�jWe4@  C0q�l���j?�0�  C/����l��j(g,�  C/	��(]��j�)   C.N�(�x�i�i%@  C-+E���i��!�  C,�B�S#n�i�k�  C,:�H�=�i��   C,�Zd��i�m@  C+�}���i���  C,�v�ȧ�ilo�  C,s�n���iT�   C,��ʴ��i=q@  C-D�j���i%��  C-X`��Wj�ir��  C-f����k�h���   C-r�w`�h�t�@  C-�A��Ț�h���  C-�J��l�h�v��  C-�¡���h���   C-����Z��h�x�@  C-�9��n�hi��  C-�7r���hRz��  C-�;hN=�h:��   C.����h#|�@  C.b��p�R�h�ր  C.���|�g�~��  C/&!�+d��g���   C/�"V}>�gŀ�@  C/��V*��g�ǀ  C/��sd�g����  C/���V�1�g�   C0	���\�gg��@  C0'}LQ�gP��  C00�}���g8���  C03������g!�   C0D(��)�g	��@  C0V�tmŘ�f�	��  C0h�ժƧ�fڊ��  C0�����P�f��   C0���)L��f���@  C0��N�z��f���  C0�w ����f|���  C1"�
}��fe�   C1+r_ޫ�fM��@  C1<�����f6��  C1O8�:CX�f���  C1a;����f�   C1+�vj��e@  C1���� �e�|�  C1
	��e��x�  C0'��X��e�u   C.��a���e��q@  C.{����ezm�  C-z_�~b��eb�i�  C-3䂰��eKf   C-	�_����e3�b@  C-
`�4Y��e^�  C-%f1%��e�Z�  C-�F���d�W   C,��kC���dՠS@  C,�R
'��d�!O�  C,�e���d��K�  C-	�����d�#H   C-�.oE�dw�D@  C-%/���
�d`%@�  C,��V�0+�dH�<�  C,�7vo���d1'9   C,���G5{�d�5@  C-J]!��d)1�  C-cks1@�c�-�  C-�Yk"�c�+*   C.0M˸��c��&@  C.�VB޼C�c�-"�  C.�Kgf�c���  C.�o��x^�cu/   C.�T��N�c]�@  C.��(\�cF1�  C.�[{�}�c.��  C/���Ս�c3   C/<)���b��@  C.��M���b�5�  C.ϵ�A��bж �  C.ߒ{�u��b�6�   C/ �k��b���@  C/:;��7��b�8��  C/�ya� k�br���  C/��l��b[:�   C0)
�3S�bC��@  C03}�oS�b,<�  C0E*d�LB�b���  C0[Sv�U��a�>�   C0w�?>��a��@  C0����A�a�@׀  C0� �NEM�a����  C0�d�ܳ*�a�B�   C0���R��a���@  C0�����L�apDȀ  C0��H��aX���  C0�:��aAF�   C1;�+(F�a)ǽ@  C1J�3���aH��  C1~6w_��`�ɵ�  C1��<���`�J�   C1ı~���`�ˮ@  C1�`�(y�`�L��  C1��i[�`�ͦ�  C2j�7���`�N�   C1��߹���`mϟ@  C/�ccL}R�`VP��  C.a�Go���`>ї�  C-PcJVkU�`'R�   C,�s����`Ӑ@  C,SX����_�   C,3�u3�_���  C,I�f:5��_��
   C,��Y�-��_c��  C-�*��_4��   C-v)���_��  C-���"�^ִ�   C-~3̆�2�^���  C-�K8�p��^x��   C-������^I�Հ  C-�Qk���^��   C-�_p B�]�ƀ  C.&�����]���   C-���b�]�·�  C-��M:�]^İ   C-���x���]/ƨ�  C-��]a��] ȡ   C."����\�ʙ�  C.z݃�y��\�̒   C.�n��Iu�\sΊ�  C/7�r4Z�\DЃ   C/I�W<�\�{�  C/o'T���[��t   C/��^@�[��l�  C/�0ƠU�[��e   C0D��Z�[Y�]�  C0B���8�[*�V   C0I�U��Z��N�  C0E!ژ�6�Z��G   C0P@�W]�Z��?�  C0`O��[�Zn�8   C0u�ٿ��Z?�0�  C0����(��Z�)   C0���OX�Y��!�  C0햤���Y��   C1�;�y:�Y���  C1:���+�YT�   C1&���
l�Y%��  C14mީ�P�X���   C1FM��Eo�X���  C1\D�����X���   C1%W��t��Xi��  C1��.�j�X:��   C1vEP��X�ր  C0�߯k��W���   C0%<�}Q��W�ǀ  C/6$���W�   C.^��c�WP��  C-����W�W!�   C-dk�[�V�	��  C-B\R�	�V��   C,�Q��ґ�V���  C,r�z	�h�Ve�   C,�I|���V6��  C+�
�.|=�V�   C+��ŵ��U�|�  C+��T	0�U�u   C+t��]W�Uzm�  C+x�c���UKf   C+��$���U^�  C+��ş���T�W   C,4�dƠ��T�!O�  C,�U�`���T�#H   C-ަ�5�T`%@�  C-�^c`i�T1'9   C-��1SS�T)1�  C-�f�^9��S�+*   C.
{��0��S�-"�  C.�hK��Su/   C.Ag��*4�SF1�  C.m�����S3   C-�'F�x�R�5�  C-�9���R�6�   C-�GA��R�8��  C-�v'�{�R[:�   C.(�Ő�#�R,<�  C.n��-��Q�>�   C.�Eg"4�Q�@׀  C/�F�0�Q�B�   C/�C ���QpDȀ  C/��|jq�QAF�   C/�)+D�QH��  C0�,���P�J�   C0"� ��.�P�L��  C0<o�>���P�N�   C0)�CNsH�PVP��  C0 -����P'R�   C0-ڰ����O�   C0H	���O��
   C0e��ڹ��O4��   C0��D�w��Nִ�   C0�5�Rе�Nx��   C0�a$���N��   C1���{��M���   C1@!w��M^İ   C1�;?���M ȡ   C08���f�L�̒   C/j�~ѭ/�LDЃ   C.�	�	��K��t   C.e���t��K��e   C-a￷�]�K*�V   C,�h
&(�J��G   C,E�i^:��Jn�8   C,.0Չ�J�)   C,��M<��I��   C,9h�W��IT�   C,o�j��t�H���   C,��bcZ��H���   C-$,��C�H:��   C-��Iq@��G���   C-�8yD��G�   C-��a���G!�   C-�1Z�0�F��   C-�,e���Fe�   C-�[4�f�F�   C.%��	�E�u   C.X��q�8�EKf   C-�!̠B��D�W   C-ܘL\��D�#H   C-�Ċ^��D1'9   C.	���n��C�+*   C.8ʗ^��Cu/   C.��P(���C3   C.���҅�B�6�   C/DC���/�B[:�   C/�ES�#a�A�>�   C/�T���A�B�   C/�T��AAF�   C02�6�-�@�J�   C0%oy(���@�N�   C0f�X�@'R�   C0j���1�?��
   C0`�z0 ��>ִ�   C0g��kzk�>��   C0y	-���=^İ   C0w�I�\��<�̒   C0Qz�l��;��t   C/5�����;*�V   C.���	��:n�8   C.Qf&ڃ�9��   C. γ{h��8���   C-u]�T���8:��   C,�����7�   C,��V�#�6��   C,4hvq>��6�   C+��w��_�5Kf   C+�"�;n��4�#H   C+��Q���3�+*   C+]�/���33   C+Lv�}���2[:�   C+T���CK�1�B�   C+^�r���0�J�   C+���P��0'R�   C,-:OS0��.ִ�   C,�k%>Mb�-^İ   C-P%~��+��t   C-ܜ1�-�*n�8   C-"'�@��(���   C-3U�j���'�   C-T�̇<��&�   C-�k��ZY�$�#H   C-���H�#3   C-z@O�P��!�B�   C-k�|u�� 'R�   C-�o�[n�^İ   C-��M��#�n�8   C-�p��v��   C.-LE�ib��#H   C.�XU|{D��B�   C.���ۓ�^İ   C/Fuz�xs��   C/X�� ����B�   C/|&h51���   C/�xB�Q���   C0Mע�        C0<�/�̓A��   C02Z����A��   C0)��j�"B�B�   C0/tb��B�   C0A9����B^İ   C0X�Qt`B�B�   C0{"�ުB�#H   C0��%2]�B�   C0�FxF�Bn�8   C0�!k,��B^İ   C1+q��NB 'R�   C1:�m�ۘB!�B�   C1IM�<e�B#3   C1T&�'��B$�#H   C0Z�\�B&�   C.E�0���B'�   C-�2xB� B(���   C,�Y	rB*n�8   C,]r��zB+��t   C,��hM�B-^İ   C,�oSuB.ִ�   C,*f��SB0'R�   C,gmdg%B0�J�   C,�_*��B1�B�   C-25�<�B2[:�   C-��}��BB33   C-���m:WB3�+*   C-���rG�B4�#H   C-�X$'�0B5Kf   C-�`A��B6�   C-� f�B6��   C.��]�jB7�   C-�!��PB8:��   C-�mwvFB8���   C-�tWW~B9��   C-�K ��B:n�8   C.*����'B;*�V   C.vt��B;��t   C.�rD��B<�̒   C/5���GB=^İ   C/�JF��JB>��   C/��ь��B>ִ�   C/Ą���B?��
   C/��}���B@'R�   C0ӾHBB@�N�   C0*,j{��B@�J�   C0G!�kBAAF�   C0'����BA�B�   C/��e�BA�>�   C.��'�BB[:�   C.�h�A�BB�6�   C-��2�O�BC3   C-0~4�!�BCu/   C-y��� BC�+*   C-D��2�BD1'9   C-"�r�BD�#H   C,�<XT��BD�W   C,�o;��PBEKf   C,n?T^BE�u   C,H^�H0cBF�   C,7o��LRBFe�   C,3���BF��   C+����5�BG!�   C+��N!4�BG�   C+�X��BG���   C+��RRBH:��   C+��u�6BH���   C+�{x�lBH���   C,RU��F�BIT�   C,���BI��   C-%��s�fBJ�)   C-4YmmڪBJn�8   C-I]�� CBJ��G   C-d���BK*�V   C-��n$9�BK��e   C-��*��BK��t   C-�Κ�[�BLDЃ   C.�^)#BL�̒   C-�^b
�wBM ȡ   C-�(��BM^İ   C.ٝ���BM���   C.5��k_�BN��   C.zEnclBNx��   C.�Au�BNִ�   C/5@�i�BO4��   C/����/BO��
   C/�R��BO�   C/� Ò�BP'R�   C/�+�8BPVP��  C/ч>���BP�N�   C0����BP�L��  C0���e�BP�J�   C0.�T��BQH��  C0(Dj�6�BQAF�   C0>R�BQpDȀ  C0Vl�\�BQ�B�   C0z${��BQ�@׀  C0�zY��#BQ�>�   C0�W�`�
BR,<�  C1��$jBR[:�   C1-!%F6�BR�8��  C1<W�%,ZBR�6�   C1PU?�BR�5�  C1d���<{BS3   C1rܘ�޿BSF1�  C1p���j�BSu/   C1E4o� �BS�-"�  C1=D��$KBS�+*   C1?��$�gBT)1�  C1F2譔%BT1'9   C1X
�IBT`%@�  C1�t3�[&BT�#H   C1���<fBT�!O�  C1�+���BT�W   C1{�/��BU^�  C0;;϶�BUKf   C/u�g& BUzm�  C.&^�SI$BU�u   C-�e�_�BU�|�  C-4۫�X�BV�   C,CM3X��BV6��  C+ف6�8�BVe�   C+�}��~BV���  C+��y�BV��   C+�}���BV�	��  C,.ݺ�&�BW!�   C,r�(��BWP��  C,�o/��BW�   C-+BW�ǀ  C-w��7��BW���   C-�Y]�5BX�ր  C-�ɏ�BX:��   C-�G��cBXi��  C.�'bb+BX���   C.�D�eBX���  C.@T�w_OBX���   C.wP"�]BY%��  C.%���kBYT�   C.	6��3�BY���  C.�w�BY��   C.:����BY��!�  C.c�.B�`BZ�)   C.��vyBZ?�0�  C.���t*BZn�8   C/d��=��BZ��?�  C/����BZ��G   C/�z�e��BZ��N�  C09��VB[*�V   C0=2jVB[Y�]�  C0b]ʒ*B[��e   C0�E�l�B[��l�  C0��dy��B[��t   C0��%E�B\�{�  C0��r�.B\DЃ   C0�c���B\sΊ�  C0���UqB\�̒   C0����N�B\�ʙ�  C0��e��B] ȡ   C0�d��qB]/ƨ�  C1)���(B]^İ   C0��,���B]�·�  C/�1&V��B]���   C.�3V��
B]�ƀ  C.])�،�B^��   C-�^����B^I�Հ  C-.92���B^x��   C,?�/o�B^���  C+�ܻ��B^ִ�   C+o��ϲ�B_��  C+J��~O~B_4��   C+C͑��/B_c��  C+p��ES�B_��
   C+�~��=�B_���  C,<��B_�   C,�����B`Ӑ@  C-��pDfB`'R�   C-"l�XB`>ї�  C-.'�KN�B`VP��  C-9��#JnB`mϟ@  C-J����B`�N�   C-n�p��#B`�ͦ�  C-��%j�B`�L��  C-ڳW�MB`�ˮ@  C-�Q�2�B`�J�   C-�m�qI�B`�ɵ�  C-���_BaH��  C-����Ba)ǽ@  C.0�@�	aBaAF�   C.����q�BaX���  C.��zV�tBapDȀ  C/FPj�gBa���@  C/W��4��Ba�B�   C/xiϴr�Ba����  C/���|Ba�@׀  C0Rz��Ba��@  C04��\��Ba�>�   C0a\�u�+Bb���  C0=Egr�Bb,<�  C0A�lߨmBbC��@  C0Oe
�E|Bb[:�   C0e]-t� Bbr���  C0���ᐃBb�8��  C0�]-JiBb���@  C0΂!�۽Bb�6�   C0��$�Bbж �  C1(&�-dZBb�5�  C19�|pr�Bb��@  C1O�	��Bc3   C1f�P*Bc.��  C/�Q��z�BcF1�  C.@���k�Bc]�@  C-.;d���Bcu/   C,V�C-��Bc���  C+��SĔBc�-"�  C+��ܤ�@Bc��&@  C+|t���Bc�+*   C+���>B�Bc�-�  C+�^}�(Bd)1�  C,N�8}$�Bd�5@  C,É��N'Bd1'9   C-&��yBdH�<�  C-3����Bd`%@�  C-9����Bdw�D@  C->��-��Bd�#H   C-Ht�w�Bd��K�  C-s?��ΖBd�!O�  C-��a�6YBdՠS@  C-�� ���Bd�W   C-�9l���Be�Z�  C-�C��hBe^�  C-����Be3�b@  C-�9^�5_BeKf   C.ìi�Beb�i�  C.qv�Ւ�Bezm�  C.��Ŝ`Be��q@  C/"�h���Be�u   C/(��xBe��x�  C/=��+:�Be�|�  C/J��K�SBe@  C/�����Bf�   C/�w�a�;Bf���  C/��+DBf6��  C/�T��BfM��@  C0D�� �Bfe�   C0X�BYxBf|���  C00cTi@Bf���  C0M����Bf���@  C0tY�=� Bf��   C0��M:Bfڊ��  C0��޿ЃBf�	��  C0��h4�YBg	��@  C1D�4�Bg!�   C1��N�Bg8���  C1*��5upBgP��  C19+T�sLBgg��@  C1[l���Bg�   C0J�p���Bg����  C.�����Bg�ǀ  C.+���5cBgŀ�@  C-��8]Bg���   C-��TqBg�~��  C,��O�Bh�ր  C,����Bh#|�@  C,�fԙ)�Bh:��   C-N���QBhRz��  C-J�#��Bhi��  C-9���}Bh�x�@  C-#ku�Bh���   C-(U�]�-Bh�v��  C-B�'Ky�Bh���  C-c���8�Bh�t�@  C-��k<�UBh���   C-��?&�Bir��  C-m�m�+�Bi%��  C-dV	�c�Bi=q@  C-y���BiT�   C-�q�@>Bilo�  C-�O�	�Bi���  C.%n"�L�Bi�m@  C.��4��Bi��   C.͏��zBi�k�  C.�r��I\Bi��!�  C.�B�Z�xBi�i%@  C/����3Bj�)   C/,�xٹ�Bj(g,�  C/_�j�WBj?�0�  C/�N��g�BjWe4@  C/��_�Bjn�8   C/���Qn?Bj�c;�  C/��=�TBj��?�  C0K���~Bj�aC@  C0)�f��Bj��G   C0S��Bj�_J�  C0�c�2yBj��N�  C0��|�Bk]R@  C0���|Bk*�V   C0�e8���BkB[Y�  C1�6D�BkY�]�  C1���<BkqYa@  C1+���(Bk��e   C1A^=Mz�Bk�Wh�  C1D���gBk��l�  C1�J�GBk�Up@  C0��,�Bk��t   C02��Dz�Bk�Sw�  C/O� $�Bl�{�  C.[ax��Bl-Q@  C-�F/�BlDЃ   C-��r�eBl\O��  C-�>6�J�BlsΊ�  C-�{��{BBl�M�@  C-@�O��Bl�̒   C-OP��Bl�K��  C,�#ʽ��Bl�ʙ�  C,ڲS�x�Bl�I�@  C,���lrBm ȡ   C,����ĐBmG��  C,]�����Bm/ƨ�  C,)Ǌ=��BmGE�@  C, ֣)�1Bm^İ   C,F��#BmvC��  C,��gfBm�·�  C,�68���Bm�A�@  C-&e�7�Bm���   C-�����Bm�?��  C-��^fHBm�ƀ  C-�ѷZ��Bn=�@  C.�	�e�Bn��   C.(��b��Bn2;��  C.A�0�BnI�Հ  C.qd���Bna9�@  C.��o^�Bnx��   C.���nvQBn�7��  C.�z�Qm�Bn���  C.���
�Bn�5�@  C.��r�y�Bnִ�   C.���<RBn�3��  C.�r�UBo��  C/C�����Bo1�@  C/���"h�Bo4��   C/��N%�BoL/��  C0}�+��Boc��  C0��Վ�Bo{.@  C00�� %ABo��
   C0I=*��Bo�,�  C0{r��Bo���  C0�ҳCxBo�*@  C0��2�ӖBo�   C0�BŖ�+Bp`  C0��>	Y�BpӐ@  C0�kI�Z