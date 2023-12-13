CDF  �   
      time       lon       lat          	   CDI       @Climate Data Interface version 2.0.3 (https://mpimet.mpg.de/cdi)   Conventions       CF-1.5     source        #PISM development no-version-control    command      / /home/m/m300792/ModelCodes/mPISM_build/install_03-10-22/bin/pismr -i pism_-065000/pism_-065000.nc -stress_balance ssa+sia -y 100 -ys -65000 -o pism_-064900/pism_-064900.nc -surface_given_file forcing/smb_-064900.nc -surface given -ocean_th_file forcing/oce_-064900.nc -ocean th -topg_file pism_-065000/pism_-065000_new_topg.nc -bed_def prescribed -uplift_file pism_-065000/pism_-065000_dbdt.nc -extra_times 10 -extra_vars thk,usurf,bheatflx,velsurf_mag,velbase_mag,tillwat,bmelt,velbase,flux_divergence,velsurf,velbar,wvelsurf,wvelbase,climatic_mass_balance_cumulative,potential_climatic_mass_balance_cumulative,tempbase,tempsurf,mask,temppabase,tempicethk_basal,topg,nonneg_flux_cumulative,grounded_basal_flux_cumulative,floating_basal_flux_cumulative,discharge_flux_cumulative,Href,taud_mag,lat,lon,calving_mask,sea_level -extra_file pism_-064900/pism_-064900_extra -extra_split -ts_times yearly -ts_file pism_-064900/pism_-064900_ts.nc -topg_to_phi 15.0,30.0,-300.0,400.0 -ssafd_ksp_rtol 1e-7 -enthalpy_temperate_conductivity_ratio 0.001 -pseudo_plastic -pseudo_plastic_q .25 -pseudo_plastic_uthreshold 70 -till_effective_fraction_overburden 0.04 -tauc_slippery_grounding_lines -sia_e 8 -ssa_e 1 -ssa_upwind_factor 0.7 -calving eigen_calving,thickness_calving -thickness_calving_threshold 100 -part_redist -part_grid -cfbc -kill_icebergs -eigen_calving_K 1e16 -verbose 2 -sliding_scale_factor_reduces_tauc 1 -th_heat_coefficient 4e-5 -th_salinity_coefficient 4e-8 -o_format netcdf3 -sediment_file till_mask.nc -no_divide_tillphi -periodicity none -ocean_th_period 1 -o_order yxz
    proj4         j+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs     history_of_appended_files        �Wed Sep 14 10:09:59 2022: Appended file outFinal2 had following "history" attribute:
Wed Sep 14 10:09:57 2022: ncap2 -O -s lon_bnds=float(lon_bnds) outFinal1 outFinal2
Wed Sep 14 10:09:57 2022: ncap2 -O -s lat_bnds=float(lat_bnds) outFinal outFinal1
Wed Sep 14 10:09:56 2022: cdo -add thkNew thkOld outFinal
Wed Sep 14 10:09:54 2022: cdo -mul -selvar,thk pism_-064900/pism_-064900.nc MaskLaurentide thkNew
     NCO       _netCDF Operators version 5.0.6 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)    history      7�Fri Sep 16 09:57:32 2022: cdo -mergetime File000001.nc File000002.nc File000003.nc File000004.nc File000005.nc File000006.nc File000007.nc File000008.nc File000009.nc File000010.nc File000011.nc File000012.nc File000013.nc File000014.nc File000015.nc File000016.nc File000017.nc File000018.nc File000019.nc File000020.nc File000021.nc File000022.nc File000023.nc File000024.nc File000025.nc File000026.nc File000027.nc File000028.nc File000029.nc File000030.nc File000031.nc File000032.nc File000033.nc File000034.nc File000035.nc File000036.nc File000037.nc File000038.nc File000039.nc File000040.nc File000041.nc File000042.nc File000043.nc File000044.nc File000045.nc File000046.nc File000047.nc File000048.nc File000049.nc File000050.nc File000051.nc File000052.nc File000053.nc File000054.nc File000055.nc File000056.nc File000057.nc File000058.nc File000059.nc File000060.nc File000061.nc File000062.nc File000063.nc File000064.nc File000065.nc File000066.nc File000067.nc File000068.nc File000069.nc File000070.nc File000071.nc File000072.nc File000073.nc File000074.nc File000075.nc File000076.nc File000077.nc File000078.nc File000079.nc File000080.nc File000081.nc File000082.nc File000083.nc File000084.nc File000085.nc File000086.nc File000087.nc File000088.nc File000089.nc File000090.nc File000091.nc File000092.nc File000093.nc File000094.nc File000095.nc File000096.nc File000097.nc File000098.nc File000099.nc File000100.nc File000101.nc File000102.nc File000103.nc File000104.nc File000105.nc File000106.nc File000107.nc File000108.nc File000109.nc File000110.nc File000111.nc File000112.nc File000113.nc File000114.nc File000115.nc File000116.nc File000117.nc File000118.nc File000119.nc File000120.nc File000121.nc File000122.nc File000123.nc File000124.nc File000125.nc File000126.nc File000127.nc File000128.nc File000129.nc File000130.nc File000131.nc File000132.nc File000133.nc File000134.nc File000135.nc File000136.nc File000137.nc File000138.nc File000139.nc File000140.nc File000141.nc File000142.nc File000143.nc File000144.nc File000145.nc File000146.nc File000147.nc File000148.nc File000149.nc File000150.nc File000151.nc File000152.nc File000153.nc File000154.nc File000155.nc File000156.nc File000157.nc File000158.nc File000159.nc File000160.nc File000161.nc File000162.nc File000163.nc File000164.nc File000165.nc File000166.nc File000167.nc File000168.nc File000169.nc File000170.nc File000171.nc File000172.nc File000173.nc File000174.nc File000175.nc File000176.nc File000177.nc File000178.nc File000179.nc File000180.nc File000181.nc File000182.nc File000183.nc File000184.nc File000185.nc File000186.nc File000187.nc File000188.nc File000189.nc File000190.nc File000191.nc File000192.nc File000193.nc File000194.nc File000195.nc File000196.nc File000197.nc File000198.nc File000199.nc File000200.nc File000201.nc File000202.nc File000203.nc File000204.nc File000205.nc File000206.nc File000207.nc File000208.nc File000209.nc File000210.nc File000211.nc File000212.nc File000213.nc File000214.nc File000215.nc File000216.nc File000217.nc File000218.nc File000219.nc File000220.nc File000221.nc File000222.nc File000223.nc File000224.nc File000225.nc File000226.nc File000227.nc File000228.nc File000229.nc File000230.nc File000231.nc File000232.nc File000233.nc File000234.nc File000235.nc File000236.nc File000237.nc File000238.nc File000239.nc File000240.nc File000241.nc File000242.nc File000243.nc File000244.nc File000245.nc File000246.nc File000247.nc File000248.nc File000249.nc File000250.nc File000251.nc File000252.nc File000253.nc File000254.nc File000255.nc File000256.nc File000257.nc File000258.nc File000259.nc File000260.nc File000261.nc File000262.nc File000263.nc File000264.nc File000265.nc File000266.nc File000267.nc File000268.nc File000269.nc File000270.nc File000271.nc File000272.nc File000273.nc File000274.nc File000275.nc File000276.nc File000277.nc File000278.nc File000279.nc File000280.nc File000281.nc File000282.nc File000283.nc File000284.nc File000285.nc File000286.nc File000287.nc File000288.nc File000289.nc File000290.nc File000291.nc File000292.nc File000293.nc File000294.nc File000295.nc File000296.nc File000297.nc File000298.nc File000299.nc File000300.nc File000301.nc File000302.nc File000303.nc File000304.nc File000305.nc File000306.nc File000307.nc File000308.nc File000309.nc File000310.nc File000311.nc File000312.nc File000313.nc File000314.nc File000315.nc File000316.nc File000317.nc File000318.nc File000319.nc File000320.nc File000321.nc File000322.nc File000323.nc File000324.nc File000325.nc File000326.nc File000327.nc File000328.nc File000329.nc File000330.nc File000331.nc File000332.nc File000333.nc File000334.nc File000335.nc File000336.nc File000337.nc File000338.nc File000339.nc File000340.nc File000341.nc File000342.nc File000343.nc File000344.nc File000345.nc File000346.nc File000347.nc File000348.nc File000349.nc File000350.nc File000351.nc File000352.nc File000353.nc File000354.nc File000355.nc File000356.nc File000357.nc File000358.nc File000359.nc File000360.nc File000361.nc File000362.nc File000363.nc File000364.nc File000365.nc File000366.nc File000367.nc File000368.nc File000369.nc File000370.nc File000371.nc File000372.nc File000373.nc File000374.nc File000375.nc File000376.nc File000377.nc File000378.nc File000379.nc File000380.nc File000381.nc File000382.nc File000383.nc File000384.nc File000385.nc File000386.nc File000387.nc File000388.nc File000389.nc File000390.nc File000391.nc File000392.nc File000393.nc File000394.nc File000395.nc File000396.nc File000397.nc File000398.nc File000399.nc File000400.nc File000401.nc File000402.nc File000403.nc File000404.nc File000405.nc File000406.nc File000407.nc File000408.nc File000409.nc File000410.nc File000411.nc File000412.nc File000413.nc File000414.nc File000415.nc File000416.nc File000417.nc File000418.nc File000419.nc File000420.nc File000421.nc File000422.nc File000423.nc File000424.nc File000425.nc File000426.nc File000427.nc File000428.nc File000429.nc File000430.nc File000431.nc File000432.nc File000433.nc File000434.nc File000435.nc File000436.nc File000437.nc File000438.nc File000439.nc File000440.nc File000441.nc File000442.nc File000443.nc File000444.nc File000445.nc File000446.nc File000447.nc File000448.nc File000449.nc File000450.nc File000451.nc File000452.nc File000453.nc File000454.nc File000455.nc File000456.nc File000457.nc File000458.nc File000459.nc File000460.nc File000461.nc File000462.nc File000463.nc File000464.nc File000465.nc File000466.nc File000467.nc File000468.nc File000469.nc File000470.nc File000471.nc File000472.nc File000473.nc File000474.nc File000475.nc File000476.nc File000477.nc File000478.nc File000479.nc File000480.nc File000481.nc File000482.nc File000483.nc File000484.nc File000485.nc File000486.nc File000487.nc File000488.nc File000489.nc File000490.nc File000491.nc File000492.nc File000493.nc File000494.nc File000495.nc File000496.nc File000497.nc File000498.nc File000499.nc File000500.nc File000501.nc File000502.nc File000503.nc File000504.nc File000505.nc File000506.nc File000507.nc File000508.nc File000509.nc File000510.nc File000511.nc File000512.nc File000513.nc File000514.nc File000515.nc File000516.nc File000517.nc File000518.nc File000519.nc File000520.nc File000521.nc File000522.nc File000523.nc File000524.nc File000525.nc File000526.nc File000527.nc File000528.nc File000529.nc File000530.nc File000531.nc File000532.nc File000533.nc File000534.nc File000535.nc File000536.nc File000537.nc File000538.nc File000539.nc File000540.nc File000541.nc File000542.nc File000543.nc File000544.nc File000545.nc File000546.nc File000547.nc File000548.nc File000549.nc File000550.nc File000551.nc File000552.nc File000553.nc File000554.nc File000555.nc File000556.nc File000557.nc File000558.nc File000559.nc File000560.nc File000561.nc File000562.nc File000563.nc File000564.nc File000565.nc File000566.nc File000567.nc File000568.nc File000569.nc File000570.nc File000571.nc File000572.nc File000573.nc File000574.nc File000575.nc File000576.nc File000577.nc File000578.nc File000579.nc File000580.nc File000581.nc File000582.nc File000583.nc File000584.nc File000585.nc File000586.nc File000587.nc File000588.nc File000589.nc File000590.nc File000591.nc File000592.nc File000593.nc File000594.nc File000595.nc File000596.nc File000597.nc File000598.nc File000599.nc File000600.nc File000601.nc File000602.nc File000603.nc File000604.nc File000605.nc File000606.nc File000607.nc File000608.nc File000609.nc File000610.nc File000611.nc File000612.nc File000613.nc File000614.nc File000615.nc File000616.nc File000617.nc File000618.nc File000619.nc File000620.nc File000621.nc File000622.nc File000623.nc File000624.nc File000625.nc File000626.nc File000627.nc File000628.nc File000629.nc File000630.nc File000631.nc File000632.nc File000633.nc File000634.nc File000635.nc File000636.nc File000637.nc File000638.nc File000639.nc File000640.nc File000641.nc File000642.nc File000643.nc File000644.nc File000645.nc File000646.nc File000647.nc File000648.nc File000649.nc File000650.nc File000651.nc File000652.nc File000653.nc File000654.nc File000655.nc File000656.nc File000657.nc File000658.nc File000659.nc File000660.nc File000661.nc File000662.nc File000663.nc File000664.nc File000665.nc File000666.nc File000667.nc File000668.nc File000669.nc File000670.nc File000671.nc File000672.nc File000673.nc File000674.nc File000675.nc File000676.nc File000677.nc File000678.nc File000679.nc File000680.nc File000681.nc File000682.nc File000683.nc File000684.nc File000685.nc File000686.nc File000687.nc File000688.nc File000689.nc File000690.nc File000691.nc File000692.nc File000693.nc File000694.nc File000695.nc File000696.nc File000697.nc File000698.nc File000699.nc File000700.nc File000701.nc File000702.nc File000703.nc File000704.nc File000705.nc File000706.nc File000707.nc File000708.nc File000709.nc File000710.nc File000711.nc File000712.nc File000713.nc File000714.nc File000715.nc File000716.nc File000717.nc File000718.nc File000719.nc File000720.nc File000721.nc File000722.nc File000723.nc File000724.nc File000725.nc File000726.nc File000727.nc File000728.nc File000729.nc File000730.nc File000731.nc File000732.nc File000733.nc File000734.nc File000735.nc File000736.nc File000737.nc File000738.nc File000739.nc File000740.nc File000741.nc File000742.nc File000743.nc File000744.nc File000745.nc File000746.nc File000747.nc File000748.nc File000749.nc File000750.nc File000751.nc File000752.nc File000753.nc File000754.nc File000755.nc File000756.nc File000757.nc File000758.nc File000759.nc File000760.nc File000761.nc File000762.nc File000763.nc File000764.nc File000765.nc File000766.nc File000767.nc File000768.nc File000769.nc File000770.nc File000771.nc File000772.nc File000773.nc File000774.nc File000775.nc File000776.nc File000777.nc File000778.nc File000779.nc File000780.nc File000781.nc File000782.nc File000783.nc File000784.nc File000785.nc File000786.nc File000787.nc File000788.nc File000789.nc File000790.nc File000791.nc File000792.nc File000793.nc File000794.nc File000795.nc File000796.nc File000797.nc File000798.nc File000799.nc File000800.nc File000801.nc File000802.nc File000803.nc File000804.nc File000805.nc File000806.nc File000807.nc File000808.nc File000809.nc File000810.nc File000811.nc File000812.nc File000813.nc File000814.nc File000815.nc File000816.nc File000817.nc File000818.nc File000819.nc File000820.nc File000821.nc File000822.nc File000823.nc File000824.nc File000825.nc File000826.nc File000827.nc File000828.nc File000829.nc File000830.nc File000831.nc File000832.nc File000833.nc File000834.nc File000835.nc File000836.nc File000837.nc File000838.nc File000839.nc File000840.nc File000841.nc File000842.nc File000843.nc File000844.nc File000845.nc File000846.nc File000847.nc File000848.nc File000849.nc File000850.nc File000851.nc File000852.nc File000853.nc File000854.nc File000855.nc File000856.nc File000857.nc File000858.nc File000859.nc File000860.nc File000861.nc File000862.nc File000863.nc File000864.nc File000865.nc File000866.nc File000867.nc File000868.nc File000869.nc File000870.nc File000871.nc File000872.nc File000873.nc File000874.nc File000875.nc File000876.nc File000877.nc File000878.nc File000879.nc File000880.nc File000881.nc File000882.nc File000883.nc File000884.nc File000885.nc File000886.nc File000887.nc File000888.nc File000889.nc File000890.nc File000891.nc File000892.nc File000893.nc File000894.nc File000895.nc File000896.nc File000897.nc File000898.nc File000899.nc File000900.nc File000901.nc File000902.nc File000903.nc File000904.nc File000905.nc File000906.nc File000907.nc File000908.nc File000909.nc File000910.nc File000911.nc File000912.nc File000913.nc File000914.nc File000915.nc File000916.nc File000917.nc File000918.nc File000919.nc File000920.nc File000921.nc File000922.nc File000923.nc File000924.nc File000925.nc File000926.nc File000927.nc File000928.nc File000929.nc File000930.nc File000931.nc File000932.nc File000933.nc File000934.nc File000935.nc File000936.nc File000937.nc File000938.nc File000939.nc File000940.nc File000941.nc File000942.nc File000943.nc File000944.nc File000945.nc File000946.nc File000947.nc File000948.nc File000949.nc File000950.nc File000951.nc File000952.nc File000953.nc File000954.nc File000955.nc File000956.nc File000957.nc File000958.nc File000959.nc File000960.nc File000961.nc File000962.nc File000963.nc File000964.nc File000965.nc File000966.nc File000967.nc File000968.nc File000969.nc File000970.nc File000971.nc File000972.nc File000973.nc File000974.nc File000975.nc File000976.nc File000977.nc File000978.nc File000979.nc File000980.nc File000981.nc File000982.nc File000983.nc File000984.nc File000985.nc File000986.nc File000987.nc File000988.nc File000989.nc File000990.nc File000991.nc File000992.nc File000993.nc File000994.nc File000995.nc File000996.nc File000997.nc File000998.nc File000999.nc File001000.nc IceVolumeHudson.nc
Fri Sep 16 09:47:59 2022: cdo -fldsum -mulc,10000 -mulc,10000 -selindexbox,279,477,230,362 -selvar,thk /work/ba0989/m300792/HE_Runs/ExperimentsComposite//HE97//pism_-064900/pism_-064900.nc Tmp/File000001.nc     CDO       @Climate Data Operators version 2.0.3 (https://mpimet.mpg.de/cdo)         time                standard_name         time   	long_name         time   units         seconds since 1-1-1    calendar      365_day    axis      T               D�   lon                standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X               D�   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y               D�   thk                       standard_name         land_ice_thickness     	long_name         land ice thickness     units         m      pism_intent       model_state             D�                �}Ȁ@�  C4����I��}���   C4������}�=   C4����}�A�@  C4ͧv�?��}��9`  C3ٶ��<�}�·�  C2�G4�tW�}�5�  C2]��W*w�}vC��  C1��,w�}j�1�  C1Q�tQ��}^İ   C0��L�(��}S.   C0ʍBT�z�}GE�@  C0�s����};�*`  C0��>���}/ƨ�  C0z���EC�}$&�  C0v�����}G��  C0����&*�}�"�  C0�hC#4�} ȡ   C0��>����|�	   C0�[����|�I�@  C1ʽ [��|݊`  C1.��W��|�ʙ�  C1H��w�|��  C1b������|�K��  C1z�U���|���  C1�1r�\�|�̒   C1������|�   C1�t20�*�|�M�@  C1��Ũy>�|�`  C1������|sΊ�  C1���Lr�|h�  C2s��9s�|\O��  C20c����|P��  C2K�����|DЃ   C2f����|9   C2��[t�|-Q@  C2�G.�x|�|!��`  C2�½��$�|�{�  C2Ɗwߚ �|
��  C2�Mڤ5�{�Sw�  C2���a���{���  C3�@��B�{��t   C3u��q�{��   C3&d� ��{�Up@  C38/zŤ�{Õ�`  C3K��Qƙ�{��l�  C3e��ږ\�{��  C3��}�p�{�Wh�  C3�CY�$��{����  C3��o��{��e   C3�;V�!�{}�   C3�����w�{qYa@  C3�-'���{e��`  C3�����{Y�]�  C4R{7�{N۠  C4�6���{B[Y�  C4,�߲/!�{6���  C4<�o?^�{*�V   C4K�����{�   C4Z�R9/�{]R@  C4i%0�0��{��`  C4z_J�B��z��N�  C4�$�X�l�z�̠  C4���+��z�_J�  C4���#g�z؟��  C4���Y<�z��G   C4Λ�Q��z� �   C3�&X*?�z�aC@  C2�Z�6(c�z���`  C29�˞���z��?�  C1�jPg���z�"��  C1=�'M(��z�c;�  C0��^wm\�zz���  C0�Y?����zn�8   C0�������zc$�   C0�������zWe4@  C0w!�ji��zK��`  C0e,F@��z?�0�  C0X0I5ȸ�z4&��  C0L��h�e�z(g,�  C0H'���c�z���  C0@�����z�)   C0Ii����z(�   C0m<R��8�y�i%@  C0�����M�y���`  C0����?�y��!�  C0��Yw7�y�*��  C0�jMh1��y�k�  C0�A�����y����  C1#Ǿ��y��   C1�ɼ��y�,�   C1=��"���y�m@  C1Uwe����y���`  C1r=�Q�g�y���  C1���j�q�yx.��  C1���$�j�ylo�  C1ü�?��y`���  C1� ��K3�yT�   C1�(�X��yI0�   C2([XKr�y=q@  C2"��@4��y1��`  C28!�j��y%��  C1�3�%���y2��  C1��ӌ��yr��  C2j��0��y�}�  C2 V���.�x���   C27B�Y@��x�4z   C2R �1��x�t�@  C2ki��xӵv`  C2���[���x���  C2������x�6r�  C2�xO���x�v��  C2�}�A�|�x��n�  C2�-wx��x���   C2�j��b\�x�8k   C3K��2��x�x�@  C3�xY} �xu�g`  C3)���@a�xi��  C3;Q��.��x^:c�  C3L)v����xRz��  C3]��d��xF�_�  C3q�!�LA�x:��   C3�>�e?�x/<\   C3��0Ž�x#|�@  C3��Cѓt�x�X`  C3˷R����x�ր  C3�)d��v�x >T�  C3�i�^�w�~��  C4T~�z�w�P�  C4�착��w���   C4)0�b���w�@M   C49��IHO�wŀ�@  C4I4E?G��w��I`  C4XWŚ��w�ǀ  C4gUC ���w�BE�  C4v�ӌc��w����  C4�6qm�y�w��A�  C4�o�-"i�w�   C4�B%̜�wsD>   C4�	c�l�wg��@  C3��G����w[�:`  C2�7�����wP��  C2����wDF6�  C1x��w8���  C1	M�| �w,�2�  C0«���E�w!�   C0��9�1z�wH/   C0u�h2X�w	��@  C0Z�X1���v��+`  C0D���u�v�	��  C05��<��v�J'�  C0+�(����vڊ��  C0$�kR�v��#�  C0ԛBʶ�v��   C0�zp/��v�L    C0�I �&�v���@  C0�la�>�v��`  C0[+�U��v���  C0#��6Ƥ�v�N�  C0@��� �v|���  C0_RK���vp��  C0x��U	�ve�   C0�����U�vYP   C0�Z[b��vM��@  C0����m��vA�`  C0�Z<�rP�v6��  C0�M��v*R	�  C1�R���v���  C1#I�y#�v��  C1A�s�~��v�   C1a��w��u�T   C1x�5��(�u@  C1���8��u���`  C1�����u�|�  C1��M�Bd�u�U��  C1�1��$�u��x�  C1�t�ч��u����  C1�Ѩ��)�u�u   C2q)�%+�u�W�   C2$����u��q@  C284w���u���`  C2N'µ�uzm�  C2j&~����unY�  C2����we�ub�i�  C2�Ӈ-U&�uV���  C2��i�e�uKf   C2ɫ����u?[�   C2߂�X݉�u3�b@  C2���C��u'��`  C3�/8���u^�  C3!�,��u]ܠ  C39����u�Z�  C3Fw�cP�t����  C3U>���t�W   C3e��T���t�_�   C3w	�|�+�tՠS@  C3��+>U��t���`  C3�\��#Z�t�!O�  C3��W��t�a͠  C3���r��t��K�  C3�H:!��t����  C3��V�V��t�#H   C4
�&<��t�c�   C4���1��tw�D@  C4/�燎��tk��`  C4@Qe'��t`%@�  C4PB��.�tTe��  C4_��j���tH�<�  C4nZp^M�t<��  C4|�ڸH�t1'9   C4�*�a���t%g�   C4��>+��t�5@  C4���/�t�`  C4Y�Z����t)1�  C3��'9�s�i��  C2$�Pb��s�-�  C1r��3���s���  C0�k�$��s�+*   C0�Z]��s�k�   C0Q�v����s��&@  C0&�m�+��s��`  C0JD[u��s�-"�  C/̅�ëY�s�m��  C/��t��a�s���  C/���,&1�s���  C/�ʽ�v��su/   C/q�d	��sio�   C/b	��3��s]�@  C/P2ݡ���sQ�`  C/@���sF1�  C/0���\r�s:q��  C/V`W#S	�s.��  C/�y�$��s"��  C/�*y�p��s3   C0���ss�   C0-M`�`�r��@  C0K@�f6�r��`  C0d���S��r�5�  C0wj`���r�u��  C0���x��rж �  C0�?{^��r��~�  C0¾�M���r�6�   C0�h۰M��r�w{   C0�`"E"�r���@  C1Ocr��r��w`  C14��&��r�8��  C1N����r~ys�  C1i"����rr���  C1�Lp�}��rf�o�  C1���j\��r[:�   C1�ʈ�lj�rO{l   C1��8&��rC��@  C1���?���r7�h`  C1�����r,<�  C2�X8V��r }d�  C2%���r���  C2-5vtC��r�`�  C2F鹀���q�>�   C2`�Ƥ���q�]   C2zQ���q��@  C2�M�,���q� Y`  C2��0�^`�q�@׀  C2����w�qU�  C2�r��#D�q����  C2�>�$��q�Q�  C2��ܣ���q�B�   C3��p���q��N   C3#"���q���@  C38�'��q|J`  C3N��Y���qpDȀ  C3aĝ���qd�F�  C3m�p4�qX���  C3��Cw~�qMB�  C3���xiq�qAF�   C3�����9�q5�?   C3�C��V�q)ǽ@  C3�xZ���q;`  C3����4�qH��  C3�DބPW�q�7�  C4�!�:D�p�ɵ�  C4 ���.s�p�
3�  C40�7&��p�J�   C4@Q���p׋0   C4N�!u���p�ˮ@  C4\2��]$�p�,`  C4j`�5�p�L��  C4xG51��p��(�  C4�PP��p�ͦ�  C2��(��i�p�$�  C1�M��%p�p�N�   C1E5���%�py�!   C0�i��a�pmϟ@  C0w�) A^�pb`  C0A�1A��pVP��  C0hB>���pJ��  C/�rv���p>ї�  C/����M�p3�  C/���u��p'R�   C/s��U�p�   C/[�^ �g�pӐ@  C/B(�$b��p`  C/'Y0����o�   C/�O]��o�*@  C.��V��o���  C/�.9��o�,�  C/LMo����o��
   C/��L��o{.@  C/ب�ݼ��oc��  C0z���oL/��  C0#���X��o4��   C0>zx�x��o1�@  C0Wm_'�f�o��  C0oJ�p�H�n�3��  C0�D�O��nִ�   C0������n�5�@  C0����%%�n���  C0�������n�7��  C0�u�O��nx��   C1a�ٯ�na9�@  C1#^�R��nI�Հ  C1>��� �n2;��  C1X�U��n��   C1q�7�2n�n=�@  C1��r���m�ƀ  C1�O���m�?��  C1�,��<��m���   C1��gNs��m�A�@  C1អc���m�·�  C1�mi�/N�mvC��  C2��\�m^İ   C2@��,&�mGE�@  C23���0$�m/ƨ�  C2M��Y��mG��  C2gQG%,-�m ȡ   C2�n�=�n�l�I�@  C2��Qt���l�ʙ�  C2�i��l�K��  C2��YWb��l�̒   C2ܘ�BQ�l�M�@  C2�z�N���lsΊ�  C3tT�;�l\O��  C3�9Q���lDЃ   C3+����l-Q@  C3<��$��l�{�  C3M�hN$�k�Sw�  C3`Fm��y�k��t   C3t������k�Up@  C3�, m��k��l�  C3��r$�k�Wh�  C3��}�j�k��e   C3��*�!�kqYa@  C3���ə��kY�]�  C3��	��kB[Y�  C4�R.���k*�V   C4���6��k]R@  C4%+Z�
p�j��N�  C3��%Ud�j�_J�  C3
�����j��G   C2V�{ �@�j�aC@  C1�Pq�_�j��?�  C1t�����j�c;�  C0���V�jn�8   C0g��x�1�jWe4@  C03|� 
��j?�0�  C0� �/��j(g,�  C/�ࢪ� �j�)   C/�c
�[�i�i%@  C/\�`�9�i��!�  C/>�L�l6�i�k�  C/B�ܱ���i��   C/;)n����i�m@  C/']E��g�i���  C/%�����ilo�  C/H����iT�   C/ ML0�i=q@  C/M����i%��  C/���5K�ir��  C/��+�n�h���   C0����h�t�@  C0&�U��n�h���  C0Ev ��h�v��  C0]�Rr��h���   C0y��5��h�x�@  C0��j�3��hi��  C0�����B�hRz��  C0ńĸ	��h:��   C0������h#|�@  C0������h�ր  C1
U�0���g�~��  C1 t�W˾�g���   C16���0�gŀ�@  C1O�5�!��g�ǀ  C1l�RU�g����  C1�r�S;��g�   C1������gg��@  C1������gP��  C1���7�?�g8���  C1�6j�fb�g!�   C2���k��g	��@  C2 d����f�	��  C2:2!�4��fڊ��  C2H+]��f��   C2[W6蝛�f���@  C2n@
�[-�f���  C2���|c��f|���  C2��k��;�fe�   C2���I�+�fM��@  C2¼���f6��  C2ۦ�\�h�f���  C2���@	�f�   C3
��]�e@  C3 �d"�e�|�  C36ܣz��e��x�  C3J��"��e�u   C3^Θ�h�e��q@  C3p�:гX�ezm�  C3�d�00�eb�i�  C3�@$���eKf   C3���	Ex�e3�b@  C3��9$
�e^�  C3���n���e�Z�  C3��h��d�W   C3��-�3?�dՠS@  C4e��� �d�!O�  C4��R��d��K�  C4,�Jd\��d�#H   C4?.�{.��dw�D@  C4MF��s��d`%@�  C3��+�G#�dH�<�  C2�aB1���d1'9   C23��F���d�5@  C1x]Հ|��d)1�  C0�v5*�h�c�-�  C0�,�X���c�+*   C0Q�}8��c��&@  C0\.Y=�c�-"�  C/�������c���  C/��i �0�cu/   C/g�,m(��c]�@  C/9Oս33�cF1�  C/#%���h�c.��  C/2����c3   C.��|���b��@  C.ܤ�e�b�5�  C.�g��L�bж �  C.��O,ͱ�b�6�   C.�|�f���b���@  C.�^L9��b�8��  C/@TKm��br���  C/X��:7p�b[:�   C/�0=���bC��@  C/� �u��b,<�  C0��MzH�b���  C0!��=���a�>�   C0@3 �x#�a��@  C0]j�s��a�@׀  C0yeЎ �a����  C0�e�2�f�a�B�   C0�?@�o�a���@  C0�:��apDȀ  C0߳6'�C�aX���  C0���vR�aAF�   C1i��c��a)ǽ@  C1!����aH��  C18j�?�`�ɵ�  C1LiX�U��`�J�   C1dNW��
�`�ˮ@  C1[x"��`�L��  C1�,;`y��`�ͦ�  C1���径�`�N�   C1��j��`mϟ@  C1��-����`VP��  C1����a��`>ї�  C2䆇P1�`'R�   C2,IY�,�`Ӑ@  C2A݂0&�_�   C2U�J����_���  C2i���q��_��
   C2|kQyS�_c��  C2�c����_4��   C2���qo�_��  C2�� \@��^ִ�   C2Ѿ�Uj��^���  C2ꔢN��^x��   C3��&��^I�Հ  C33C#���^��   C30����]�ƀ  C3EE�>���]���   C3Y�����]�·�  C3m{��M�]^İ   C3�c #D��]/ƨ�  C3��'�#.�] ȡ   C3�g��\�ʙ�  C3�D,�Q�\�̒   C3������\sΊ�  C3�����\DЃ   C3���U*�\�{�  C3�M+'8,�[��t   C45O���[��l�  C4%�>Qy�[��e   C42?�i�[Y�]�  C3r�WD�f�[*�V   C2��zD�Z��N�  C1�2VR��Z��G   C1{�� �Z��?�  C0�/��3T�Zn�8   C0L�m_���Z?�0�  C0���8�Z�)   C/ӪЦ,;�Y��!�  C/��8}���Y��   C/X,�>�Y���  C/(��U&�YT�   C.�z�l���Y%��  C.���N�X���   C.�cs��X���  C.��&�d:�X���   C.Գ�H�r�Xi��  C.ǊF�^��X:��   C.�*�2��X�ր  C.��s�&�W���   C.��bu�W�ǀ  C/:u$���W�   C/s�'�)�WP��  C/�hI����W!�   C/�.LΖ�V�	��  C/�ײ�Z�V��   C0�<B�V���  C08��˃��Ve�   C0W*�*���V6��  C0t���Y��V�   C0�z{���U�|�  C0��d����U�u   C0�h���q�Uzm�  C0�P�(P�UKf   C0���3W��U^�  C1��*���T�W   C1*E�ڢ�T�!O�  C1?����{�T�#H   C1Sa���#�T`%@�  C1h���J*�T1'9   C1}�i�T��T)1�  C1��?�XS�S�+*   C1��1K��S�-"�  C1�S��_�Su/   C1��͍P�SF1�  C2�ˬ�Z�S3   C2�(h�|�R�5�  C22��:\�R�6�   C2I�tX��R�8��  C2`3yu
.�R[:�   C2u��ŉ��R,<�  C2�� K�j�Q�>�   C2�(��ݬ�Q�@׀  C2��&1�n�Q�B�   C2�T�ۆ��QpDȀ  C2�Z6��QAF�   C2���e՛�QH��  C3|����P�J�   C3�\+_��P�L��  C36�����P�N�   C3MT��(�PVP��  C3cG�IV��P'R�   C3v�չz�O�   C3��/���O��
   C3��a����O4��   C3�@E�1�Nִ�   C3��Cz�Nx��   C3�4�yK�N��   C3�8H{�"�M���   C3���oD�M^İ   C4��*/c�M ȡ   C4q�z���L�̒   C4(>�}��LDЃ   C4;����T�K��t   C4G���eF�K��e   C3m�ß��K*�V   C2{��Y�J��G   C1�b�I��Jn�8   C0���7K��J�)   C0���Y4�I��   C0;�G.�[�IT�   C0�9�x�H���   C/��1`<�H���   C/��];|]�H:��   C/Q�:Z��G���   C//<#΂_�G�   C/QdI��G!�   C.�p�����F��   C.��	���Fe�   C.�+�l�}�F�   C.�qe�Ý�E�u   C.�VJ�!�EKf   C.����qs�D�W   C.�~���w�D�#H   C.��lӬ��D1'9   C.��@g���C�+*   C/:ë�7^�Cu/   C/n<'���C3   C/�=@��K�B�6�   C/�Y�U?�B[:�   C0_��"�A�>�   C0#P0��A�B�   C0=�)+���AAF�   C0[W�2��@�J�   C0y��@�N�   C0����!��@'R�   C0�өI0�?��
   C0�S�{~�>ִ�   C0��9�C�>��   C0��٤�=^İ   C1(� (W�<�̒   C12rC`�j�;��t   C1J+�/		�;*�V   C1WJ#�M�:n�8   C1j���9��   C1������8���   C1�@L����8:��   C1��b��7�   C1̰4ņ��6��   C1����9�6�   C2 K[�.��5Kf   C2ČF�x�4�#H   C20A\����3�+*   C2G�A�L�33   C2]  �^�2[:�   C2r�v���1�B�   C2�N�ɠ_�0�J�   C2��W�^��0'R�   C2�05�l��.ִ�   C2�h���-^İ   C2�se�c�+��t   C2�~^W=�*n�8   C3 �1�r~�(���   C3�
�$��'�   C31eu����&�   C3Hh����$�#H   C3^G"����#3   C3sP�pI��!�B�   C3�np�[� 'R�   C3��o`��^İ   C3��%~��n�8   C3�_��Ȓ��   C3�L�	�3��#H   C3�]�t�*��B�   C3��XX��^İ   C3��_�Q��   C4d�����B�   C4%K�):����   C45D4��i���   C3䰑0ϵ        C2�Q���A��   C1�ҹ~A��   C16W"2-�B�B�   C0�]��
bB�   C0V��B^İ   C0�Z˦�B�B�   C/�l���/B�#H   C/h��)�B�   C/5�?տBn�8   C/#�TI�B^İ   C.��y�ہB 'R�   C.�$d�+�B!�B�   C.�їь$B#3   C.�'a�B$�#H   C.�8V�|B&�   C.����B'�   C/IE|��B(���   C/`{UB*n�8   C/0�gb�B+��t   C/q���s<B-^İ   C/��7��B.ִ�   C/�^�wB0'R�   C0�Q0B0�J�   C0$���zB1�B�   C0:ԡmn�B2[:�   C0Tq��bB33   C0r�ņ��B3�+*   C0�d\hB4�#H   C0��ސ�B5Kf   C0ĵ�ʹ1B6�   C0������B6��   C0�����B7�   C1�hA�B8:��   C1.��qU�B8���   C1H�4�B9��   C1bi���6B:n�8   C1v唜իB;*�V   C1�Vn���B;��t   C1�@g��B<�̒   C1�%�_B=^İ   C1��].�B>��   C1�~R�3B>ִ�   C1��6���B?��
   C2�.pO�B@'R�   C2+�����B@�N�   C2D]����B@�J�   C2\T�>{BAAF�   C2r��v��BA�B�   C2�����BA�>�   C2��o�SBB[:�   C2���UZ}BB�6�   C2�m��dtBC3   C2�~���
BCu/   C2�d�}�IBC�+*   C2���`�BD1'9   C3G�B��BD�#H   C3"tbi�BD�W   C3:���h�BEKf   C3R�i�MBE�u   C3i��QBF�   C3a��_BFe�   C3�Y���BF��   C3�6~%2BG!�   C3����BG�   C3̫�L�BG���   C3��u��BH:��   C3�6r�hBH���   C3�T����BH���   C4
dp!BIT�   C4#"2��BI��   C4,\�zTGBJ�)   C4=��3E�BJn�8   C4M�{�#BJ��G   C3��˹��BK*�V   C2�욘OBK��e   C1�ů�c*BK��t   C1.��rBLDЃ   C0� 0 �BL�̒   C0F9B_s�BM ȡ   C0�b0y�BM^İ   C/ĄKBP3BM���   C/�Q���sBN��   C/g~���FBNx��   C/F�C�}�BNִ�   C/4���GBO4��   C/(T0ڢBO��
   C/&�;�BO�   C/"�ؙxkBP'R�   C/%!?d;BPVP��  C/0� ���BP�N�   C/?� ��BP�L��  C/\�8�k�BP�J�   C/�����@BQH��  C/�FD�qBQAF�   C0y���BQpDȀ  C0%(~��DBQ�B�   C0=J�4#BQ�@׀  C0Ut�XBQ�>�   C0k!�˕BR,<�  C0�����YBR[:�   C0�*�-'BR�8��  C0���q��BR�6�   C0Ѐ�^.BR�5�  C0�B�B�VBS3   C1]��]BSF1�  C1#th�ԸBSu/   C1>@FG@�BS�-"�  C1W���A9BS�+*   C1r��H\BT)1�  C1�FsO�BT1'9   C1��"�ZBT`%@�  C1��E\uBT�#H   C1μ��8BT�!O�  C1��{}BT�W   C1�0�BU^�  C2m��dBUKf   C2��=�BUzm�  C27'0�eBU�u   C2Q1=��BU�|�  C2jJ����BV�   C2��:�BV6��  C2�{�sMBVe�   C2�FV[-:BV���  C2�2�Ƌ�BV��   C2�.���BV�	��  C2�O@���BW!�   C3xB�J�BWP��  C3��0��BW�   C3'��BW�ǀ  C39�J�
BW���   C3K���BX�ր  C3_�,W��BX:��   C3v��x�RBXi��  C3�����BX���   C3�pU�1BX���  C3��u�IBX���   C3�G!ɆdBY%��  C3�O#(��BYT�   C3�n��ڞBY���  C4����BY��   C4�T��+BY��!�  C4#���B�BZ�)   C41���BZ?�0�  C4=��
�BZn�8   C4K��2�BZ��?�  C4lg��BZ��G   C2�E=vBZ��N�  C2����B[*�V   C1ULB1�B[Y�]�  C0ȳ��yB[��e   C0Z��!1B[��l�  C09@�B[��t   C/�F]�D�B\�{�  C/~WB��B\DЃ   C/P.�`ʗB\sΊ�  C/1�Q'�B\�̒   C/!i��nB\�ʙ�  C/D��^B] ȡ   C/�L;.vB]/ƨ�  C/-xc�yB]^İ   C/$�H�B]�·�  C/<\d��B]���   C/�� g�+B]�ƀ  C/�:n��B^��   C0	0��=)B^I�Հ  C0(�{oB^x��   C0D!"t�B^���  C0_Ǆ�1|B^ִ�   C0z<B���B_��  C0�Hu� �B_4��   C0��ٯ�MB_c��  C0�)�[eWB_��
   C0��d�>B_���  C0�Wfi�B_�   C1	��s�B`Ӑ@  C1 ja ��B`'R�   C19���}*B`>ї�  C1U�ȝ#�B`VP��  C1rG�8�B`mϟ@  C1����?B`�N�   C1���L�B`�ͦ�  C1��~7�B`�L��  C1�*�4�B`�ˮ@  C1�qwxB`�J�   C2JP;B`�ɵ�  C2���eBaH��  C23p�Ba)ǽ@  C2E���r{BaAF�   C2Y"�ΐBaX���  C2k.͙�BapDȀ  C2}�t�~Ba���@  C2�"�E�<Ba�B�   C2�@`�r{Ba����  C2���TBa�@׀  C2��NݟBa��@  C2�]�V�Ba�>�   C3����UBb���  C3 p2w��Bb,<�  C35�P���BbC��@  C3H��%��Bb[:�   C3[���|Bbr���  C3nnϥ�Bb�8��  C3.����Bb���@  C3���r�Bb�6�   C3�|
���Bbж �  C3�Z�{ZBb�5�  C3����lBb��@  C3۬nk�pBc3   C3�d�fBc.��  C4��@�BcF1�  C4T�q5Bc]�@  C4%�X[��Bcu/   C3d]oΆ�Bc���  C2ww.*W�Bc�-"�  C1��WbzaBc��&@  C1�&��}Bc�+*   C0����Bc�-�  C0K3��t�Bd)1�  C0fGC��Bd�5@  C/����yBd1'9   C/sŉ�9BdH�<�  C/7~ےw~Bd`%@�  C/��|��Bdw�D@  C.�:/E�Bd�#H   C.�Q�TBd��K�  C.�ˍmBd�!O�  C.�T��"BdՠS@  C.�b�>!Bd�W   C.�PV��Be�Z�  C.~���c�Be^�  C.�Y����Be3�b@  C.�p�u%4BeKf   C.�D_ABeb�i�  C/4ry�Bezm�  C/^��ͯBe��q@  C/�U��i�Be�u   C/��ɓ��Be��x�  C0��>�Be�|�  C0 �Z"��Be@  C0>�W�;1Bf�   C0\qTK�Bf���  C0waRnU�Bf6��  C0��k�#BfM��@  C0��"wZ
Bfe�   C0�3|���Bf|���  C0ށiߞBf���  C0�D���Bf���@  C1��d��Bf��   C1!��a�Bfڊ��  C18��C��Bf�	��  C1SOq9'�Bg	��@  C1pA��Bg!�   C1��L$��Bg8���  C1�� ���BgP��  C1�ŗpn�Bgg��@  C1�#�[ABg�   C1�al@8_Bg����  C2�u�\Bg�ǀ  C2V��=Bgŀ�@  C25��>�Bg���   C2Ko�]�_Bg�~��  C2`7�ᶣBh�ր  C2tM#>�1Bh#|�@  C2���C1Bh:��   C2�!+o�BhRz��  C2���#�Bhi��  C2���5�Bh�x�@  C2�a�ܹBh���   C2�AN�Bh�v��  C3	\ l��Bh���  C3 ���`Bh�t�@  C35���#Bh���   C3K5Qg��Bir��  C3_i)�.Bi%��  C3r<�M��Bi=q@  C3��I�q�BiT�   C3�f�C#Bilo�  C3�:���Bi���  C3�NPS�Bi�m@  C3�i��c5Bi��   C3�̐���Bi�k�  C3�傴~�Bi��!�  C4 %9]�Bi�i%@  C4j��?Bj�)   C4+�y�a�Bj(g,�  C4==�.̙Bj?�0�  C4H�{�BjWe4@  C3ee�Ȅ�Bjn�8   C2v�=��Bj�c;�  C1��'�Bj��?�  C0�ʩh�Bj�aC@  C0�$��TKBj��G   C09I���Bj�_J�  C0_\�+�Bj��N�  C/��4�ekBk]R@  C/z��ƕBk*�V   C/>(
�f�BkB[Y�  C/2#ָ*BkY�]�  C.�#�� BkqYa@  C.۪�x�Bk��e   C.��ӕe�Bk�Wh�  C.��eU�(Bk��l�  C.h>�M�3Bk�Up@  C.KN^��Bk��t   C..u�4ܓBk�Sw�  C.(G(!:CBl�{�  C.6�� ��Bl-Q@  C.Il��t�BlDЃ   C.^Y�[Bl\O��  C.�l��BlsΊ�  C.���ȘBl�M�@  C/H �Bl�̒   C/@�<��:Bl�K��  C/~��E��Bl�ʙ�  C/�����Bl�I�@  C/�ƝB&'Bm ȡ   C0_��qBmG��  C01��_Bm/ƨ�  C0M�~���BmGE�@  C0i���H�Bm^İ   C0}ڸ�ycBmvC��  C0��4�iJBm�·�  C0�I��/�Bm�A�@  C0��ِmBm���   C0�B6dh�Bm�?��  C0����Bm�ƀ  C1��e'!Bn=�@  C1"r��w2Bn��   C1>���tBn2;��  C1Z��Q�~BnI�Հ  C1x����Bna9�@  C1��H�Bnx��   C1�w4���Bn�7��  C1�|R�%XBn���  C1����Bn�5�@  C1� ���gBnִ�   C2Q&�]�Bn�3��  C2,�+��Bo��  C2+���m%Bo1�@  C2=؄��&Bo4��   C2QzHʥYBoL/��  C2g�gw<`Boc��  C2�����Bo{.@  C2�Jؐ�Bo��
   C2�;s�Bo�,�  C2�,�B�Bo���  C2�E.�Bo�*@  C2��W�p�Bo�   C3�VBp`  C3"dց�\BpӐ@  C35�"��