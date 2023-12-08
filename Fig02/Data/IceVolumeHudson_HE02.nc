CDF  �   
      time       lon       lat          	   CDI       @Climate Data Interface version 2.0.3 (https://mpimet.mpg.de/cdi)   Conventions       CF-1.5     source        $PISM development v0.7.3-165-g7780f00   command      & /work/ba0989/m300019/mPISM-c-staggered/install/bin/pismr -i pism_-065000/pism_-065000.nc -stress_balance ssa+sia -y 100 -ys -65000 -o pism_-064900/pism_-064900.nc -surface_given_file forcing/smb_-064900.nc -surface given -ocean_th_file forcing/oce_-064900.nc -ocean th -topg_file pism_-065000/pism_-065000_new_topg.nc -bed_def prescribed -uplift_file pism_-065000/pism_-065000_dbdt.nc -extra_times 10 -extra_vars thk,usurf,bheatflx,velsurf_mag,velbase_mag,tillwat,bmelt,velbase,flux_divergence,velsurf,velbar,wvelsurf,wvelbase,climatic_mass_balance_cumulative,potential_climatic_mass_balance_cumulative,tempbase,tempsurf,mask,temppabase,tempicethk_basal,topg,nonneg_flux_cumulative,grounded_basal_flux_cumulative,floating_basal_flux_cumulative,discharge_flux_cumulative,Href,taud_mag,lat,lon,calving_mask,sea_level -extra_file pism_-064900/pism_-064900_extra -extra_split -ts_times yearly -ts_file pism_-064900/pism_-064900_ts.nc -topg_to_phi 15.0,30.0,-300.0,400.0 -ssafd_ksp_rtol 1e-7 -enthalpy_temperate_conductivity_ratio 0.001 -pseudo_plastic -pseudo_plastic_q .25 -pseudo_plastic_uthreshold 70 -till_effective_fraction_overburden 0.04 -tauc_slippery_grounding_lines -sia_e 8 -ssa_e 1 -ssa_upwind_factor 0.7 -calving eigen_calving,thickness_calving -thickness_calving_threshold 100 -part_redist -part_grid -cfbc -kill_icebergs -eigen_calving_K 1e16 -verbose 2 -sliding_scale_factor_reduces_tauc 1 -th_heat_coefficient 4e-5 -th_salinity_coefficient 4e-8 -o_format netcdf3 -sediment_file till_mask.nc -no_divide_tillphi -periodicity none -ocean_th_period 1 -o_order yxz
     proj4         j+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs     history_of_appended_files        �Wed Aug 25 10:36:52 2021: Appended file outFinal2 had following "history" attribute:
Wed Aug 25 10:36:49 2021: ncap2 -O -s lon_bnds=float(lon_bnds) outFinal1 outFinal2
Wed Aug 25 10:36:49 2021: ncap2 -O -s lat_bnds=float(lat_bnds) outFinal outFinal1
Wed Aug 25 10:36:48 2021: cdo -add thkNew thkOld outFinal
Wed Aug 25 10:36:47 2021: cdo -mul -selvar,thk pism_-064900/pism_-064900.nc MaskLaurentide thkNew
     NCO       20210825   history      D>Tue Sep 13 14:49:59 2022: cdo -mergetime File000001.nc File000002.nc File000003.nc File000004.nc File000005.nc File000006.nc File000007.nc File000008.nc File000009.nc File000010.nc File000011.nc File000012.nc File000013.nc File000014.nc File000015.nc File000016.nc File000017.nc File000018.nc File000019.nc File000020.nc File000021.nc File000022.nc File000023.nc File000024.nc File000025.nc File000026.nc File000027.nc File000028.nc File000029.nc File000030.nc File000031.nc File000032.nc File000033.nc File000034.nc File000035.nc File000036.nc File000037.nc File000038.nc File000039.nc File000040.nc File000041.nc File000042.nc File000043.nc File000044.nc File000045.nc File000046.nc File000047.nc File000048.nc File000049.nc File000050.nc File000051.nc File000052.nc File000053.nc File000054.nc File000055.nc File000056.nc File000057.nc File000058.nc File000059.nc File000060.nc File000061.nc File000062.nc File000063.nc File000064.nc File000065.nc File000066.nc File000067.nc File000068.nc File000069.nc File000070.nc File000071.nc File000072.nc File000073.nc File000074.nc File000075.nc File000076.nc File000077.nc File000078.nc File000079.nc File000080.nc File000081.nc File000082.nc File000083.nc File000084.nc File000085.nc File000086.nc File000087.nc File000088.nc File000089.nc File000090.nc File000091.nc File000092.nc File000093.nc File000094.nc File000095.nc File000096.nc File000097.nc File000098.nc File000099.nc File000100.nc File000101.nc File000102.nc File000103.nc File000104.nc File000105.nc File000106.nc File000107.nc File000108.nc File000109.nc File000110.nc File000111.nc File000112.nc File000113.nc File000114.nc File000115.nc File000116.nc File000117.nc File000118.nc File000119.nc File000120.nc File000121.nc File000122.nc File000123.nc File000124.nc File000125.nc File000126.nc File000127.nc File000128.nc File000129.nc File000130.nc File000131.nc File000132.nc File000133.nc File000134.nc File000135.nc File000136.nc File000137.nc File000138.nc File000139.nc File000140.nc File000141.nc File000142.nc File000143.nc File000144.nc File000145.nc File000146.nc File000147.nc File000148.nc File000149.nc File000150.nc File000151.nc File000152.nc File000153.nc File000154.nc File000155.nc File000156.nc File000157.nc File000158.nc File000159.nc File000160.nc File000161.nc File000162.nc File000163.nc File000164.nc File000165.nc File000166.nc File000167.nc File000168.nc File000169.nc File000170.nc File000171.nc File000172.nc File000173.nc File000174.nc File000175.nc File000176.nc File000177.nc File000178.nc File000179.nc File000180.nc File000181.nc File000182.nc File000183.nc File000184.nc File000185.nc File000186.nc File000187.nc File000188.nc File000189.nc File000190.nc File000191.nc File000192.nc File000193.nc File000194.nc File000195.nc File000196.nc File000197.nc File000198.nc File000199.nc File000200.nc File000201.nc File000202.nc File000203.nc File000204.nc File000205.nc File000206.nc File000207.nc File000208.nc File000209.nc File000210.nc File000211.nc File000212.nc File000213.nc File000214.nc File000215.nc File000216.nc File000217.nc File000218.nc File000219.nc File000220.nc File000221.nc File000222.nc File000223.nc File000224.nc File000225.nc File000226.nc File000227.nc File000228.nc File000229.nc File000230.nc File000231.nc File000232.nc File000233.nc File000234.nc File000235.nc File000236.nc File000237.nc File000238.nc File000239.nc File000240.nc File000241.nc File000242.nc File000243.nc File000244.nc File000245.nc File000246.nc File000247.nc File000248.nc File000249.nc File000250.nc File000251.nc File000252.nc File000253.nc File000254.nc File000255.nc File000256.nc File000257.nc File000258.nc File000259.nc File000260.nc File000261.nc File000262.nc File000263.nc File000264.nc File000265.nc File000266.nc File000267.nc File000268.nc File000269.nc File000270.nc File000271.nc File000272.nc File000273.nc File000274.nc File000275.nc File000276.nc File000277.nc File000278.nc File000279.nc File000280.nc File000281.nc File000282.nc File000283.nc File000284.nc File000285.nc File000286.nc File000287.nc File000288.nc File000289.nc File000290.nc File000291.nc File000292.nc File000293.nc File000294.nc File000295.nc File000296.nc File000297.nc File000298.nc File000299.nc File000300.nc File000301.nc File000302.nc File000303.nc File000304.nc File000305.nc File000306.nc File000307.nc File000308.nc File000309.nc File000310.nc File000311.nc File000312.nc File000313.nc File000314.nc File000315.nc File000316.nc File000317.nc File000318.nc File000319.nc File000320.nc File000321.nc File000322.nc File000323.nc File000324.nc File000325.nc File000326.nc File000327.nc File000328.nc File000329.nc File000330.nc File000331.nc File000332.nc File000333.nc File000334.nc File000335.nc File000336.nc File000337.nc File000338.nc File000339.nc File000340.nc File000341.nc File000342.nc File000343.nc File000344.nc File000345.nc File000346.nc File000347.nc File000348.nc File000349.nc File000350.nc File000351.nc File000352.nc File000353.nc File000354.nc File000355.nc File000356.nc File000357.nc File000358.nc File000359.nc File000360.nc File000361.nc File000362.nc File000363.nc File000364.nc File000365.nc File000366.nc File000367.nc File000368.nc File000369.nc File000370.nc File000371.nc File000372.nc File000373.nc File000374.nc File000375.nc File000376.nc File000377.nc File000378.nc File000379.nc File000380.nc File000381.nc File000382.nc File000383.nc File000384.nc File000385.nc File000386.nc File000387.nc File000388.nc File000389.nc File000390.nc File000391.nc File000392.nc File000393.nc File000394.nc File000395.nc File000396.nc File000397.nc File000398.nc File000399.nc File000400.nc File000401.nc File000402.nc File000403.nc File000404.nc File000405.nc File000406.nc File000407.nc File000408.nc File000409.nc File000410.nc File000411.nc File000412.nc File000413.nc File000414.nc File000415.nc File000416.nc File000417.nc File000418.nc File000419.nc File000420.nc File000421.nc File000422.nc File000423.nc File000424.nc File000425.nc File000426.nc File000427.nc File000428.nc File000429.nc File000430.nc File000431.nc File000432.nc File000433.nc File000434.nc File000435.nc File000436.nc File000437.nc File000438.nc File000439.nc File000440.nc File000441.nc File000442.nc File000443.nc File000444.nc File000445.nc File000446.nc File000447.nc File000448.nc File000449.nc File000450.nc File000451.nc File000452.nc File000453.nc File000454.nc File000455.nc File000456.nc File000457.nc File000458.nc File000459.nc File000460.nc File000461.nc File000462.nc File000463.nc File000464.nc File000465.nc File000466.nc File000467.nc File000468.nc File000469.nc File000470.nc File000471.nc File000472.nc File000473.nc File000474.nc File000475.nc File000476.nc File000477.nc File000478.nc File000479.nc File000480.nc File000481.nc File000482.nc File000483.nc File000484.nc File000485.nc File000486.nc File000487.nc File000488.nc File000489.nc File000490.nc File000491.nc File000492.nc File000493.nc File000494.nc File000495.nc File000496.nc File000497.nc File000498.nc File000499.nc File000500.nc File000501.nc File000502.nc File000503.nc File000504.nc File000505.nc File000506.nc File000507.nc File000508.nc File000509.nc File000510.nc File000511.nc File000512.nc File000513.nc File000514.nc File000515.nc File000516.nc File000517.nc File000518.nc File000519.nc File000520.nc File000521.nc File000522.nc File000523.nc File000524.nc File000525.nc File000526.nc File000527.nc File000528.nc File000529.nc File000530.nc File000531.nc File000532.nc File000533.nc File000534.nc File000535.nc File000536.nc File000537.nc File000538.nc File000539.nc File000540.nc File000541.nc File000542.nc File000543.nc File000544.nc File000545.nc File000546.nc File000547.nc File000548.nc File000549.nc File000550.nc File000551.nc File000552.nc File000553.nc File000554.nc File000555.nc File000556.nc File000557.nc File000558.nc File000559.nc File000560.nc File000561.nc File000562.nc File000563.nc File000564.nc File000565.nc File000566.nc File000567.nc File000568.nc File000569.nc File000570.nc File000571.nc File000572.nc File000573.nc File000574.nc File000575.nc File000576.nc File000577.nc File000578.nc File000579.nc File000580.nc File000581.nc File000582.nc File000583.nc File000584.nc File000585.nc File000586.nc File000587.nc File000588.nc File000589.nc File000590.nc File000591.nc File000592.nc File000593.nc File000594.nc File000595.nc File000596.nc File000597.nc File000598.nc File000599.nc File000600.nc File000601.nc File000602.nc File000603.nc File000604.nc File000605.nc File000606.nc File000607.nc File000608.nc File000609.nc File000610.nc File000611.nc File000612.nc File000613.nc File000614.nc File000615.nc File000616.nc File000617.nc File000618.nc File000619.nc File000620.nc File000621.nc File000622.nc File000623.nc File000624.nc File000625.nc File000626.nc File000627.nc File000628.nc File000629.nc File000630.nc File000631.nc File000632.nc File000633.nc File000634.nc File000635.nc File000636.nc File000637.nc File000638.nc File000639.nc File000640.nc File000641.nc File000642.nc File000643.nc File000644.nc File000645.nc File000646.nc File000647.nc File000648.nc File000649.nc File000650.nc File000651.nc File000652.nc File000653.nc File000654.nc File000655.nc File000656.nc File000657.nc File000658.nc File000659.nc File000660.nc File000661.nc File000662.nc File000663.nc File000664.nc File000665.nc File000666.nc File000667.nc File000668.nc File000669.nc File000670.nc File000671.nc File000672.nc File000673.nc File000674.nc File000675.nc File000676.nc File000677.nc File000678.nc File000679.nc File000680.nc File000681.nc File000682.nc File000683.nc File000684.nc File000685.nc File000686.nc File000687.nc File000688.nc File000689.nc File000690.nc File000691.nc File000692.nc File000693.nc File000694.nc File000695.nc File000696.nc File000697.nc File000698.nc File000699.nc File000700.nc File000701.nc File000702.nc File000703.nc File000704.nc File000705.nc File000706.nc File000707.nc File000708.nc File000709.nc File000710.nc File000711.nc File000712.nc File000713.nc File000714.nc File000715.nc File000716.nc File000717.nc File000718.nc File000719.nc File000720.nc File000721.nc File000722.nc File000723.nc File000724.nc File000725.nc File000726.nc File000727.nc File000728.nc File000729.nc File000730.nc File000731.nc File000732.nc File000733.nc File000734.nc File000735.nc File000736.nc File000737.nc File000738.nc File000739.nc File000740.nc File000741.nc File000742.nc File000743.nc File000744.nc File000745.nc File000746.nc File000747.nc File000748.nc File000749.nc File000750.nc File000751.nc File000752.nc File000753.nc File000754.nc File000755.nc File000756.nc File000757.nc File000758.nc File000759.nc File000760.nc File000761.nc File000762.nc File000763.nc File000764.nc File000765.nc File000766.nc File000767.nc File000768.nc File000769.nc File000770.nc File000771.nc File000772.nc File000773.nc File000774.nc File000775.nc File000776.nc File000777.nc File000778.nc File000779.nc File000780.nc File000781.nc File000782.nc File000783.nc File000784.nc File000785.nc File000786.nc File000787.nc File000788.nc File000789.nc File000790.nc File000791.nc File000792.nc File000793.nc File000794.nc File000795.nc File000796.nc File000797.nc File000798.nc File000799.nc File000800.nc File000801.nc File000802.nc File000803.nc File000804.nc File000805.nc File000806.nc File000807.nc File000808.nc File000809.nc File000810.nc File000811.nc File000812.nc File000813.nc File000814.nc File000815.nc File000816.nc File000817.nc File000818.nc File000819.nc File000820.nc File000821.nc File000822.nc File000823.nc File000824.nc File000825.nc File000826.nc File000827.nc File000828.nc File000829.nc File000830.nc File000831.nc File000832.nc File000833.nc File000834.nc File000835.nc File000836.nc File000837.nc File000838.nc File000839.nc File000840.nc File000841.nc File000842.nc File000843.nc File000844.nc File000845.nc File000846.nc File000847.nc File000848.nc File000849.nc File000850.nc File000851.nc File000852.nc File000853.nc File000854.nc File000855.nc File000856.nc File000857.nc File000858.nc File000859.nc File000860.nc File000861.nc File000862.nc File000863.nc File000864.nc File000865.nc File000866.nc File000867.nc File000868.nc File000869.nc File000870.nc File000871.nc File000872.nc File000873.nc File000874.nc File000875.nc File000876.nc File000877.nc File000878.nc File000879.nc File000880.nc File000881.nc File000882.nc File000883.nc File000884.nc File000885.nc File000886.nc File000887.nc File000888.nc File000889.nc File000890.nc File000891.nc File000892.nc File000893.nc File000894.nc File000895.nc File000896.nc File000897.nc File000898.nc File000899.nc File000900.nc File000901.nc File000902.nc File000903.nc File000904.nc File000905.nc File000906.nc File000907.nc File000908.nc File000909.nc File000910.nc File000911.nc File000912.nc File000913.nc File000914.nc File000915.nc File000916.nc File000917.nc File000918.nc File000919.nc File000920.nc File000921.nc File000922.nc File000923.nc File000924.nc File000925.nc File000926.nc File000927.nc File000928.nc File000929.nc File000930.nc File000931.nc File000932.nc File000933.nc File000934.nc File000935.nc File000936.nc File000937.nc File000938.nc File000939.nc File000940.nc File000941.nc File000942.nc File000943.nc File000944.nc File000945.nc File000946.nc File000947.nc File000948.nc File000949.nc File000950.nc File000951.nc File000952.nc File000953.nc File000954.nc File000955.nc File000956.nc File000957.nc File000958.nc File000959.nc File000960.nc File000961.nc File000962.nc File000963.nc File000964.nc File000965.nc File000966.nc File000967.nc File000968.nc File000969.nc File000970.nc File000971.nc File000972.nc File000973.nc File000974.nc File000975.nc File000976.nc File000977.nc File000978.nc File000979.nc File000980.nc File000981.nc File000982.nc File000983.nc File000984.nc File000985.nc File000986.nc File000987.nc File000988.nc File000989.nc File000990.nc File000991.nc File000992.nc File000993.nc File000994.nc File000995.nc File000996.nc File000997.nc File000998.nc File000999.nc File001000.nc File001001.nc File001002.nc File001003.nc File001004.nc File001005.nc File001006.nc File001007.nc File001008.nc File001009.nc File001010.nc File001011.nc File001012.nc File001013.nc File001014.nc File001015.nc File001016.nc File001017.nc File001018.nc File001019.nc File001020.nc File001021.nc File001022.nc File001023.nc File001024.nc File001025.nc File001026.nc File001027.nc File001028.nc File001029.nc File001030.nc File001031.nc File001032.nc File001033.nc File001034.nc File001035.nc File001036.nc File001037.nc File001038.nc File001039.nc File001040.nc File001041.nc File001042.nc File001043.nc File001044.nc File001045.nc File001046.nc File001047.nc File001048.nc File001049.nc File001050.nc File001051.nc File001052.nc File001053.nc File001054.nc File001055.nc File001056.nc File001057.nc File001058.nc File001059.nc File001060.nc File001061.nc File001062.nc File001063.nc File001064.nc File001065.nc File001066.nc File001067.nc File001068.nc File001069.nc File001070.nc File001071.nc File001072.nc File001073.nc File001074.nc File001075.nc File001076.nc File001077.nc File001078.nc File001079.nc File001080.nc File001081.nc File001082.nc File001083.nc File001084.nc File001085.nc File001086.nc File001087.nc File001088.nc File001089.nc File001090.nc File001091.nc File001092.nc File001093.nc File001094.nc File001095.nc File001096.nc File001097.nc File001098.nc File001099.nc File001100.nc File001101.nc File001102.nc File001103.nc File001104.nc File001105.nc File001106.nc File001107.nc File001108.nc File001109.nc File001110.nc File001111.nc File001112.nc File001113.nc File001114.nc File001115.nc File001116.nc File001117.nc File001118.nc File001119.nc File001120.nc File001121.nc File001122.nc File001123.nc File001124.nc File001125.nc File001126.nc File001127.nc File001128.nc File001129.nc File001130.nc File001131.nc File001132.nc File001133.nc File001134.nc File001135.nc File001136.nc File001137.nc File001138.nc File001139.nc File001140.nc File001141.nc File001142.nc File001143.nc File001144.nc File001145.nc File001146.nc File001147.nc File001148.nc File001149.nc File001150.nc File001151.nc File001152.nc File001153.nc File001154.nc File001155.nc File001156.nc File001157.nc File001158.nc File001159.nc File001160.nc File001161.nc File001162.nc File001163.nc File001164.nc File001165.nc File001166.nc File001167.nc File001168.nc File001169.nc File001170.nc File001171.nc File001172.nc File001173.nc File001174.nc File001175.nc File001176.nc File001177.nc File001178.nc File001179.nc File001180.nc File001181.nc File001182.nc File001183.nc File001184.nc File001185.nc File001186.nc File001187.nc File001188.nc File001189.nc File001190.nc File001191.nc File001192.nc File001193.nc File001194.nc File001195.nc File001196.nc File001197.nc File001198.nc File001199.nc File001200.nc File001201.nc File001202.nc File001203.nc File001204.nc File001205.nc File001206.nc File001207.nc File001208.nc File001209.nc File001210.nc File001211.nc File001212.nc File001213.nc File001214.nc File001215.nc File001216.nc File001217.nc File001218.nc File001219.nc File001220.nc File001221.nc File001222.nc File001223.nc File001224.nc File001225.nc File001226.nc File001227.nc File001228.nc File001229.nc File001230.nc IceVolumeHudson.nc
Tue Sep 13 14:39:42 2022: cdo -fldsum -mulc,10000 -mulc,10000 -selindexbox,279,477,230,362 -selvar,thk /scratch/m/m300792/HEDownload/HE02/HE02//pism_-064900/pism_-064900.nc Tmp/File000001.nc     CDO       @Climate Data Operators version 2.0.3 (https://mpimet.mpg.de/cdo)         time                standard_name         time   	long_name         time   units         seconds since 1-1-1    calendar      365_day    axis      T               Q   lon                standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X               P�   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y               P�   thk                       standard_name         land_ice_thickness     	long_name         land ice thickness     units         m      pism_intent       model_state             Q                �}Ȁ@�  C4��m���}���   C4�˶D��}�=   C4��D#���}�A�@  C4�����}��9`  C4�U҂�N�}�·�  C4Te�W�	�}�5�  C3A�MeB��}vC��  C2�U|��v�}j�1�  C2 I:w��}^İ   C1��m$��}S.   C1�r�+��}GE�@  C0�܊u��};�*`  C0�%�T���}/ƨ�  C0��3i�}$&�  C0�Sh��}G��  C0rig�3�}�"�  C0ohf�k�} ȡ   C0��Ip��|�	   C0������|�I�@  C0�)Jg��|݊`  C0��o���|�ʙ�  C0�B��S��|��  C1ؓ�4J�|�K��  C1.&<M���|���  C1Gwd����|�̒   C1`$�Q3��|�   C1x�Y�mh�|�M�@  C1�n�G��|�`  C1��t�d��|sΊ�  C1�kZ����|h�  C1ڐ����|\O��  C1�)#�|P��  C2
W)!��|DЃ   C2"*J,:@�|9   C29͠
��|-Q@  C2Q�}��|!��`  C2hm�U��|�{�  C2
��Z��|
��  C2�1�|j��{�Sw�  C2�X8Q���{���  C2�����{��t   C2��8\��{��   C2�Kg����{�Up@  C3�4�o��{Õ�`  C3��3lo�{��l�  C3+�S����{��  C3@N��\v�{�Wh�  C3T�T�M��{����  C3i"T��$�{��e   C3}!d��{}�   C3�U���{qYa@  C3����{e��`  C3����	��{Y�]�  C3�j㯢�{N۠  C3ޖ��-�{B[Y�  C3�Q�ԧ��{6���  C4����{*�V   C4�-hH�{�   C4'����{]R@  C49��CP�{��`  C4K�U o�z��N�  C4\{��=��z�̠  C4nG+P��z�_J�  C4;z����z؟��  C4����z��G   C4���;��z� �   C4��r��z�aC@  C4��zX���z���`  C4�ҤQ��z��?�  C4�+�@X��z�"��  C4<�	�C��z�c;�  C3#�F�-�zz���  C2k��6��zn�8   C1��ş><�zc$�   C1P�y1���zWe4@  C1���`�zK��`  C0Ѓk�rh�z?�0�  C0�0�;A�z4&��  C0����xJ�z(g,�  C0x��1�v�z���  C0iz�.g�z�)   C0\���,��z(�   C0_
$]�&�y�i%@  C0|)�@�]�y���`  C0�"CR�8�y��!�  C0�ԯ��y�*��  C0ٹ�n��y�k�  C0�k�N�X�y����  C1C��:�y��   C1*E�^��y�,�   C1C֚�0�y�m@  C1]7����y���`  C1u��+��y���  C1�O|3a�yx.��  C1�Gl�4�ylo�  C1�_iO�y`���  C1قs��-�yT�   C1��ȉl��yI0�   C2л��y=q@  C2 :T2���y1��`  C26�����y%��  C2M�"
t��y2��  C2d�N,]��yr��  C2{S�"��y�}�  C2�r�!��x���   C2�.�h��x�4z   C2����0��x�t�@  C2ҏ�v�h�xӵv`  C2�~�F��x���  C2��]x#��x�6r�  C3T�_[��x�v��  C3&+=���x��n�  C39�V�1��x���   C3M:���_�x�8k   C3`�B�>�x�x�@  C3t�6���xu�g`  C3�����xi��  C3�Q��`�x^:c�  C3�7wu&��xRz��  C3��	L`��xF�_�  C3�����x:��   C3�E�e���x/<\   C3�iqbJ�x#|�@  C4~�T�q�x�X`  C4 ��$���x�ր  C42�$1��x >T�  C4E$�@n3�w�~��  C4WEM���w�P�  C4fįi��w���   C4z��բ��w�@M   C4�n2S��wŀ�@  C4��"��w��I`  C4���bU��w�ǀ  C4���h.�w�BE�  C4�*��#�w����  C2���4�K�w��A�  C2@4Y~���w�   C1��xC���wsD>   C1#���61�wg��@  C0Վ��w[�:`  C0���S���wP��  C0n��Ð��wDF6�  C0Jdz`��w8���  C00�����w,�2�  C0O��W��w!�   C0�x�S�wH/   C/믌MU��w	��@  C/�`:�� �v��+`  C/壱�=��v�	��  C/���Q��v�J'�  C/�N�(���vڊ��  C/�Rcv���v��#�  C/���@�v��   C/ ,���v�L    C0(
�# �v���@  C0 �g���v��`  C0=�lŪ��v���  C0V�4RE�v�N�  C0p=Q���v|���  C0��NR���vp��  C0����ve�   C0�����<�vYP   C0����@��vM��@  C0�`~+�vA�`  C1*����v6��  C1$0�*/��v*R	�  C1=�Ӿ��v���  C1U�~��a�v��  C1n=g�y��v�   C1����;��u�T   C1���A�g�u@  C1��
v��u���`  C1�Al�"��u�|�  C1��)k��u�U��  C1����4��u��x�  C2"2[��u����  C2-i�b[�u�u   C2D����u�W�   C2[���.��u��q@  C2r�|�5�u���`  C2�_&�\��uzm�  C2�@S�Ue�unY�  C2�"l���ub�i�  C2�y�-;��uV���  C2ߢ��P�uKf   C2���?9��u?[�   C3	�a��?�u3�b@  C3$��U��u'��`  C34rc[�o�u^�  C3H����u]ܠ  C3]������u�Z�  C3q����.�t����  C3��rؽ8�t�W   C3���@@��t�_�   C3�K�Q���tՠS@  C3�����t���`  C3������t�!O�  C3�P�P���t�a͠  C3�I<�*�t��K�  C4�-�3��t����  C4 $K�3��t�#H   C41��	���t�c�   C4D�����tw�D@  C4U}$63Z�tk��`  C4f�>�C3�t`%@�  C4v����O�tTe��  C4�)19B��tH�<�  C4�d&.��t<��  C4�F2����t1'9   C4��i-���t%g�   C3�g�Gs��t�5@  C2��=���t�`  C1��$�3�t)1�  C1Nh�rN�s�i��  C0��lq��s�-�  C0��S�E��s���  C0����)��s�+*   C0n/�[��s�k�   C0R:4ٕ�s��&@  C0=�>�aa�s��`  C0,��Ȝ �s�-"�  C0#T1��c�s�m��  C0 �����s���  C0)�&;B�s���  C0J9i@���su/   C0j��?Ȝ�sio�   C0���{��s]�@  C0�Wy�!��sQ�`  C0���x���sF1�  C0���y��s:q��  C0���D9V�s.��  C1���j�s"��  C1(ɷ�N}�s3   C1B�����ss�   C1[5]쬌�r��@  C1t�_8��r��`  C1����ZU�r�5�  C1���L��r�u��  C1�9�����rж �  C1֓=����r��~�  C1�ŤZ�r�6�   C24,L�r�w{   C2j���r���@  C25���Qt�r��w`  C2L�:�� �r�8��  C2c@��k*�r~ys�  C2y����rr���  C2�v9���rf�o�  C2���f��r[:�   C2���?��rO{l   C2�<h"�rC��@  C2������r7�h`  C2�����r,<�  C3�j$u��r }d�  C3$�1��f�r���  C38���fR�r�`�  C3K�Ltؑ�q�>�   C3a;ntEq�q�]   C3z��9��q��@  C3�ePR%��q� Y`  C3�K�(��q�@׀  C3� )�p�qU�  C3���^���q����  C3ӱ���w�q�Q�  C3�
�^�|�q�B�   C3�a	�5��q��N   C4
��k\�q���@  C4vJ�q|J`  C40k��qpDȀ  C4A�(nל�qd�F�  C4R��}���qX���  C4d"�SV��qMB�  C4u#U��o�qAF�   C4��� ���q5�?   C4�����q)ǽ@  C4��|��X�q;`  C3�m��b��qH��  C2�i�~Dc�q�7�  C1��z+2��p�ɵ�  C1]�vɡx�p�
3�  C0�:R7��p�J�   C0�`�,��p׋0   C0R7H��[�p�ˮ@  C0'i�k��p�,`  C0�� u�p�L��  C/�l�hP �p��(�  C/���j:��p�ͦ�  C/��r^E�p�$�  C/qM$+��p�N�   C/S�&���py�!   C/@����pmϟ@  C/+Mn���pb`  C/��o��pVP��  C/�I�a��pJ��  C/QL�G��p>ї�  C/���]Y�p3�  C/٠�g Z�p'R�   C0����3�p�   C0#Y)�M��pӐ@  C0=��t���p`  C0Z�� ��o�   C0tw���o�*@  C0�p%�>E�o���  C0�\Ak��o�,�  C0��h��o��
   C0�[ԑL��o{.@  C0�ǿM�oc��  C1"���V�oL/��  C1,T!S� �o4��   C1>Od�X�o1�@  C1V�2�+�o��  C1o�|���n�3��  C1�<��s�nִ�   C1�S[���n�5�@  C1�0��*b�n���  C1��\�E��n�7��  C1��l٣r�nx��   C1�Y��d�na9�@  C2��I�'�nI�Հ  C2-݂����n2;��  C2D�l��<�n��   C2[�}8��n=�@  C2r��2���m�ƀ  C2�"e���m�?��  C2��߲���m���   C2��{5X��m�A�@  C2�>ǜ#3�m�·�  C2ބLX�5�mvC��  C2�w�4�m^İ   C3QCi�&�mGE�@  C3������m/ƨ�  C32X)���mG��  C3Fٺ����m ȡ   C3Z���߷�l�I�@  C3o��[
s�l�ʙ�  C3�r�4a}�l�K��  C3���2#�l�̒   C3���7	�l�M�@  C3����B�lsΊ�  C3Р"\��l\O��  C3��>U���lDЃ   C3������l-Q@  C4	h
\�l�{�  C4�l�C�k�Sw�  C4,2���k��t   C4=4I
/�k�Up@  C4LM|�u�k��l�  C4V�R����k�Wh�  C3�i�b5�k��e   C2��q>[E�kqYa@  C1��|?��kY�]�  C1K�4�=��kB[Y�  C0˅�,FF�k*�V   C0xr@B���k]R@  C0<&)&��j��N�  C0m&�yE�j�_J�  C/��e��j��G   C/��
��j�aC@  C/x#��P�j��?�  C/S�>3�j�c;�  C/J��v �jn�8   C/I1fL�3�jWe4@  C/3W٥TF�j?�0�  C/"��P}�j(g,�  C/K�4��j�)   C.��]6�d�i�i%@  C/�|4��i��!�  C/W�ᷪ��i�k�  C/��*�Y�i��   C/��]�!�i�m@  C0�&��i���  C0$��y��ilo�  C0@e��aD�iT�   C0[�IF���i=q@  C0v�V���i%��  C0�=��&��ir��  C0��@�x�h���   C0����h�t�@  C0���}$�h���  C0������h�v��  C1v�����h���   C1-�Ȅ�h�x�@  C1FC�~U�hi��  C1_8��.L�hRz��  C1w�5+�i�h:��   C1����k�h#|�@  C1�+y���h�ր  C1�1����g�~��  C1��5��-�g���   C1��tS�gŀ�@  C25Yؐ�g�ǀ  C2T��jI�g����  C26I�m�g�   C2L��`�gg��@  C2bXm��gP��  C2w�D��g8���  C2�B�q���g!�   C2��_����g	��@  C2��.�3Q�f�	��  C2���.�fڊ��  C2�E��c�f��   C2�<p���f���@  C3��l#x�f���  C3 ��`��f|���  C352�1�:�fe�   C3I�Y���fM��@  C3]��u��f6��  C3qr ?�f���  C3�2k�{�f�   C3�+ͦ���e@  C3���:%�e�|�  C3�S�6P��e��x�  C3�$7���e�u   C3�J��^��e��q@  C3�*e���ezm�  C4���Q��eb�i�  C4yZӮ.�eKf   C4&��5�p�e3�b@  C48[��T��e^�  C4F�=p��e�Z�  C3������d�W   C2������dՠS@  C2-�э��d�!O�  C1o���p��d��K�  C0�!uF�v�d�#H   C0��P����dw�D@  C0JD8IGy�d`%@�  C07a���dH�<�  C/��D��d1'9   C/�*���u�d�5@  C/q�/��+�d)1�  C/PV��H�c�-�  C/=e�_W�c�+*   C/(��[�c��&@  C/9�]�c�-"�  C.�dk)�c���  C.���=�u�cu/   C.̯����c]�@  C.��}i&��cF1�  C.���m�c.��  C/6��o��c3   C/{q|���b��@  C/��g|���b�5�  C/�:�X�bж �  C0LOf��b�6�   C0+��"7.�b���@  C0E�W#��b�8��  C0`��(4�br���  C0z�����b[:�   C0���*��bC��@  C0����b,<�  C0�H��p��b���  C0���!V�a�>�   C0�2ZX2��a��@  C1`��y�a�@׀  C1-(�B:��a����  C1F;�˲��a�B�   C1^ndll�a���@  C1v�����apDȀ  C1�7_�	X�aX���  C1�k���aAF�   C1�D��F?�a)ǽ@  C1�u;d �aH��  C1�ڸ��T�`�ɵ�  C2Уfs��`�J�   C2WJa�`�ˮ@  C25�/����`�L��  C2Lp
�A�`�ͦ�  C2a��d��`�N�   C2w³��!�`mϟ@  C2�������`VP��  C2�j&����`>ї�  C2���X!�`'R�   C2Ε`~���`Ӑ@  C2�GaGk�_�   C2���L}�_���  C30�e��_��
   C3#�}��!�_c��  C38[[����_4��   C3My_'Lq�_��  C3b=����^ִ�   C3u�݁�F�^���  C3�.�y��^x��   C3��˅%�^I�Հ  C3�}�T���^��   C3Ĳ��_�]�ƀ  C3׌a�_��]���   C3�] y	�]�·�  C3����L��]^İ   C4�uV��]/ƨ�  C4 2���_�] ȡ   C4.�#���\�ʙ�  C4;�?�2G�\�̒   C3�t�#j�\sΊ�  C2�3
���\DЃ   C1��!4��\�{�  C19d���[��t   C0���h��[��l�  C0am�p��[��e   C0$��[y�[Y�]�  C/�I8Q͏�[*�V   C/�����Z��N�  C/n]�p��Z��G   C/J�����Z��?�  C//�yI�n�Zn�8   C/������Z?�0�  C.�4����Z�)   C.�l:ݘ��Y��!�  C.�
=��^�Y��   C.«����Y���  C.�`�N��YT�   C.�h�$��Y%��  C/1����X���   C/c��� �X���  C/��f�/��X���   C/�`h����Xi��  C0U�����X:��   C0 �����X�ր  C0<Z,���W���   C0WD�0;��W�ǀ  C0q�9\/�W�   C0�i��Nv�WP��  C0����ͳ�W!�   C0��@]�\�V�	��  C0�藅�h�V��   C0��Z�V��V���  C1
c"��Ve�   C1'��ѡW�V6��  C1A	R���V�   C1ZOøv,�U�|�  C1s3J���U�u   C1��a��Uzm�  C1��}���UKf   C1���հ�U^�  C1�q��[�T�W   C1��y����T�!O�  C2p��8�T�#H   C2����T`%@�  C234�����T1'9   C2J"��QH�T)1�  C2`)��	�S�+*   C2u�����S�-"�  C2��*��Su/   C2�l��~v�SF1�  C2��e����S3   C2��o��R�5�  C2����,i�R�6�   C2���u�[�R�8��  C3
���1�R[:�   C3z���i�R,<�  C32�k�\��Q�>�   C3G.�J	~�Q�@׀  C3[}Nзe�Q�B�   C3n���&C�QpDȀ  C3����f�QAF�   C3�/Ұ)��QH��  C3�șk��P�J�   C3�������P�L��  C3�?���#�P�N�   C3�î� ��PVP��  C3�ml�,B�P'R�   C4 -���O�   C4��\��O��
   C3~G�6���O4��   C2����Z�Nִ�   C1�@��J�Nx��   C16i&v�A�N��   C0���&�M���   C0]���#�M^İ   C0"z�Z�M ȡ   C/�5KU�A�L�̒   C/�@���LDЃ   C/n�-�m�K��t   C/= �ig�K��e   C/�,V�K*�V   C/0��1;�J��G   C/�"�oI�Jn�8   C.�S���J�)   C.݅�����I��   C.�����IT�   C.�d�%��H���   C.�,^?�H���   C.�7�����H:��   C.��y�H?�G���   C.���P�G�   C/!=�>���G!�   C/`J��W�F��   C/��\�*��Fe�   C/�@4�!�F�   C0�ojJ��E�u   C0������EKf   C08{O�7��D�W   C0Ri�M^��D�#H   C0m��1��D1'9   C0�D#~�C�+*   C0���h�Cu/   C0��w��\�C3   C0�lc��B�6�   C0�vGS/��B[:�   C1�i����A�>�   C1z
T�a�A�B�   C18G��M�AAF�   C1Q?�	��@�J�   C1i��&��@�N�   C1����3��@'R�   C1�$��#v�?��
   C1�YӃ
�>ִ�   C1�TH��A�>��   C1��Z��=^İ   C1��J)�F�<�̒   C2� KY�;��t   C2(`��!��;*�V   C2?��ߨ�:n�8   C2VCX.��9��   C2ly*j�8���   C2��W��{�8:��   C2�{-|�4�7�   C2�"�+���6��   C2±�0[��6�   C2���)�
�5Kf   C2�<�ū��4�#H   C3�{����3�+*   C3�ҳ��33   C3,饵���2[:�   C3@븴��1�B�   C3V U����0�J�   C3j]<~�B�0'R�   C3}ϸT�>�.ִ�   C3�D�/�-^İ   C3�^w�d,�+��t   C3�w)U"9�*n�8   C3ˋ����(���   C3�3�$�{�'�   C3�ZP{Q��&�   C4���p�$�#H   C4c�}
�#3   C4!�/ǚ��!�B�   C4'�\}�� 'R�   C34��v#��^İ   C2\��(Gt�n�8   C1�?r�1Z��   C0�M1h���#H   C0��� ����B�   C08M(�^İ   C0��TP��   C/�t�Ff���B�   C/z3���J���   C/?��ʲ����   C/�қ�A        C/}��m�A��   C.�ɇ-&A��   C.��t��OB�B�   C.�ݥ#B+B�   C.����B^İ   C.����B�B�   C.��23o�B�#H   C.�Nl�o�B�   C.ُ��GBn�8   C/!�?�B^İ   C/d���cB 'R�   C/���SB!�B�   C/�|,BB#3   C0P��� B$�#H   C0!�hF�B&�   C0=IR~7�B'�   C0X`�>O�B(���   C0s�/c�B*n�8   C0��-_�B+��t   C0���DB-^İ   C0��cI�B.ִ�   C0�:�:B0'R�   C0��<��B0�J�   C1d�D��B1�B�   C1(��ږB2[:�   C1A�ߛ��B33   C1Z��UVfB3�+*   C1s3����B4�#H   C1���8�KB5Kf   C1��~Aq�B6�   C1���b�uB6��   C1Ӫ!��B7�   C1�5.}�SB8:��   C2���a�B8���   C2ċ�;B9��   C20Ɲ��B:n�8   C2Gۙg� B;*�V   C2^l=D�'B;��t   C2tT�I�dB<�̒   C2���T��B=^İ   C2�����B>��   C2��#ԨFB>ִ�   C2��m&��B?��
   C2���U��B@'R�   C2�}�tB@�N�   C3���$KB@�J�   C3n��BAAF�   C30�DG[SBA�B�   C3E[�4BA�>�   C3YE�_qBB[:�   C3lb0׆BB�6�   C3��-��dBC3   C3�T�7�BCu/   C3���/��BC�+*   C3�m�=��BD1'9   C3�P0l{�BD�#H   C3�S��jeBD�W   C3�~�p�BEKf   C3m�i��BE�u   C2�Bɿ�;BF�   C1�I;��BFe�   C16�5��BF��   C0���a6BG!�   C0V�4��mBG�   C0�)jecBG���   C/͚8'�	BH:��   C/����sBH���   C/Q���YBH���   C/ �5ZbBIT�   C.�e���BI��   C.�5(F�BJ�)   C.֊�o/�BJn�8   C.���fBJ��G   C.�n�4�@BK*�V   C.�Ï5�BK��e   C.�F�#BK��t   C.�ϙU13BLDЃ   C.�J�Q_BL�̒   C.�5���BM ȡ   C/C3�e�BM^İ   C/[�$H�/BM���   C/�"h�a�BN��   C/���4"BNx��   C/�&��Z�BNִ�   C0�
+2�BO4��   C01���BO��
   C0LZ���BO�   C0d��:;BP'R�   C0~i���BPVP��  C0�K�NBP�N�   C0�+ؔ��BP�L��  C0���u��BP�J�   C0�HmmR^BQH��  C0�+�xBQAF�   C1F�j�_BQpDȀ  C10����>BQ�B�   C1Jk_9BQ�@׀  C1b���BQ�>�   C1y��qnBR,<�  C1�/�[��BR[:�   C1�x})Q�BR�8��  C1�>�Ѫ�BR�6�   C1�3^�l�BR�5�  C1�'�ɯ�BS3   C2A���BSF1�  C2"��޺JBSu/   C29�M_ɬBS�-"�  C2Q
���BS�+*   C2hH�?MBT)1�  C23:m�BT1'9   C2�ꞷM`BT`%@�  C2�b]�1>BT�#H   C2¬�T�BT�!O�  C2���^��BT�W   C2��H:�BU^�  C3Y1BUKf   C3�9�ܫBUzm�  C3,�z���BU�u   C3A����BU�|�  C3U_��BV�   C3i��	BV6��  C3|kQ��4BVe�   C3��@��BV���  C3�ﾼpBV��   C3�r�J(BV�	��  C3�ݪ�3BW!�   C3�r�/�BWP��  C3�N�BW�   C4 ��)BW�ǀ  C4(����BW���   C4$R�<_BX�ր  C45)3��1BX:��   C4B)��b+BXi��  C3�d�ɳ BX���   C2�x:�c�BX���  C2K�z�BX���   C1V��hG�BY%��  C0ĺ��|BYT�   C0j��s0�BY���  C0)Ur�BY��   C/���CBY��!�  C/��"��BZ�)   C/C�j���BZ?�0�  C/�KTZBZn�8   C.�d.8��BZ��?�  C.�^P�K}BZ��G   C.�y�V�2BZ��N�  C.�c��,�B[*�V   C.ܟ�!�B[Y�]�  C.㥮�"B[��e   C.�Q(.�9B[��l�  C/|,�=B[��t   C/2��1B\�{�  C/{�(r�B\DЃ   C/�p��B\sΊ�  C/�W�0��B\�̒   C0!^�B\�ʙ�  C0.M���B] ȡ   C0L�|�B]/ƨ�  C0iec���B]^İ   C0{$]M�B]�·�  C0����9B]���   C0��W�XB]�ƀ  C0�]�'?B^��   C0�~�'kB^I�Հ  C0�-���B^x��   C1V߲k B^���  C1*��*�B^ִ�   C1B��.DB_��  C1[n/�WLB_4��   C1s���9�B_c��  C1��GߦB_��
   C1��\��B_���  C1�p��ûB_�   C1��,�!B`Ӑ@  C1�%Q��B`'R�   C2K�I�%B`>ї�  C2eu�"WB`VP��  C2/@YpXB`mϟ@  C2E��I�lB`�N�   C2\`
8�`B`�ͦ�  C2r�_���B`�L��  C2����W�B`�ˮ@  C2��quB`�J�   C2�y����B`�ɵ�  C2�"�V2BaH��  C2��ޤ�`Ba)ǽ@  C2�R �BaAF�   C3
+�rnmBaX���  C3`"�ҠBapDȀ  C33�s��Ba���@  C3H�U/kBa�B�   C3]c���kBa����  C3p%�xj�Ba�@׀  C3��ˡl]Ba��@  C3��Z�iBa�>�   C3����Bb���  C3�W���Bb,<�  C3�#���BbC��@  C3�v�s��Bb[:�   C3򥘡%TBbr���  C4�1�$QBb�8��  C4�����Bb���@  C4$D<�Bb�6�   C4��ۊ#Bbж �  C2���\�Bb�5�  C2,{̒BBb��@  C1v�)G�Bc3   C0��$5nBc.��  C0��l���BcF1�  C0.��e{@Bc]�@  C/��>�MNBcu/   C/�-�)2�Bc���  C/PȒ��Bc�-"�  C/"�U��Bc��&@  C/��[�Bc�+*   C.��v���Bc�-�  C.��dh�WBd)1�  C.�z�(�Bd�5@  C.����Bd1'9   C/�O1�BdH�<�  C/*�>W�LBd`%@�  C/pM�d�Bdw�D@  C/��hZA�Bd�#H   C/�Ct��Bd��K�  C0/' 2�Bd�!O�  C0)�Iy��BdՠS@  C0G���Bd�W   C0d��p�sBe�Z�  C0v�}���Be^�  C0�<���Be3�b@  C0��7��eBeKf   C0�9^��Beb�i�  C0�'e��Bezm�  C0�
��Be��q@  C1��!Be�u   C1$��+s�Be��x�  C1=����Be�|�  C1V)t�{Be@  C1n��qBf�   C1�*mj0Bf���  C1���f<Bf6��  C1���%BfM��@  C1�����Bfe�   C1����kBf|���  C1���j�Bf���  C2�W��Bf���@  C2.NG�_�Bf��   C2Eh)k�Bfڊ��  C2\�����Bf�	��  C2s���Bg	��@  C2��\i�Bg!�   C2�m�% 3Bg8���  C2��\>~BgP��  C2̽�=Bgg��@  C2�V�4Bg�   C2�'Y�u'Bg����  C3Ӳ^
�Bg�ǀ  C3"n$�ZBgŀ�@  C37Ch�?�Bg���   C3K3�>E�Bg�~��  C3`#�:�Bh�ր  C3sx"�#jBh#|�@  C3�!z&�Bh:��   C3��6�BhRz��  C3�R n�Bhi��  C3�7�M�vBh�x�@  C3�&�EgBh���   C3�R��۸Bh�v��  C3����7�Bh���  C4B���Bh�t�@  C4���әBh���   C4$�'�Bir��  C4$b
Э�Bi%��  C3��4Bi=q@  C2?=���@BiT�   C1�l�*3Bilo�  C0�x�]hcBi���  C0��B��Bi�m@  C02���a/Bi��   C/�p�A�]Bi�k�  C/��d���Bi��!�  C/cꅾb�Bi�i%@  C/9��Bj�)   C/(�� �Bj(g,�  C/�V�Bj?�0�  C/�EF��BjWe4@  C/~PPhBjn�8   C/�֩uBj�c;�  C/f��"�Bj��?�  C/(���IBj�aC@  C/B����Bj��G   C/���Bj�_J�  C/�:4���Bj��N�  C0�,��Bk]R@  C0ʩA��Bk*�V   C09��֫NBkB[Y�  C0U_auBkY�]�  C0o���W�BkqYa@  C0���ʮBk��e   C0�!{b��Bk�Wh�  C0��a�4�Bk��l�  C0���L�Bk�Up@  C0�&�V$�Bk��t   C1
����SBk�Sw�  C1"ЫOE�Bl�{�  C1;(�pJ1Bl-Q@  C1R�t�BlDЃ   C1k��'��Bl\O��  C1����x�BlsΊ�  C1��w�.Bl�M�@  C1�Q�CoPBl�̒   C1��/Bl�K��  C1�`%�XBl�ʙ�  C1���j��Bl�I�@  C2�2+�LBm ȡ   C2'���BmG��  C2>��c�[Bm/ƨ�  C2U���sBmGE�@  C2k�r�_�Bm^İ   C2�;+�#BmvC��  C2���syBm�·�  C2����ZBm�A�@  C2ë7�ѰBm���   C2������Bm�?��  C2�:[�N�Bm�ƀ  C3X��^Bn=�@  C3t���QBn��   C3-e��GBn2;��  C3A�4�pBnI�Հ  C3U���ABna9�@  C3jq��{4Bnx��   C3~�8���Bn�7��  C3����E)Bn���  C3�-�!#Bn�5�@  C3�� 7��Bnִ�   C3���&ԛBn�3��  C3�;�ٚ�Bo��  C3��a���Bo1�@  C4gz���Bo4��   C4�cg�BoL/��  C4&�`��Boc��  C43��ŭ%Bo{.@  C3��0֗JBo��
   C2��e�Bo�,�  C1���-��Bo���  C1N�(�+mBo�*@  C0�n�Tx�Bo�   C0e��>�Bp`  C09D[�~BpӐ@  C/§>ƴIBp�   C/��`wZBp'R�   C/Q����Bp3�  C/-�]TblBp>ї�  C/�i ��BpJ��  C/X)��RBpVP��  C/
��Bpb`  C/
؍g1�Bpmϟ@  C/M�)��Bpy�!   C/�N���Bp�N�   C/?�ι�Bp�$�  C/�!d�y Bp�ͦ�  C/���EBp��(�  C/�s�V�Bp�L��  C0_3�X�Bp�,`  C02/$�WBp�ˮ@  C0I����Bp׋0   C0f �Bp�J�   C0��_Bp�
3�  C0����Bp�ɵ�  C0�Q&B�Bq�7�  C0�8�s�nBqH��  C0�-�xۚBq;`  C0��V.Bq)ǽ@  C1x�p0_Bq5�?   C11���5%BqAF�   C1J��WABqMB�  C1cb��BqX���  C1{0I���Bqd�F�  C1��Q�BqpDȀ  C1�/]��Bq|J`  C1��	�S�Bq���@  C1��+�Bq��N   C1�w⍭Bq�B�   C2����Bq�Q�  C2k[�9Bq����  C26�1��BqU�  C2M��˭�Bq�@׀  C2dY#t�aBq� Y`  C2zƔ�T�Bq��@  C2��_"@Bq�]   C2���q��Bq�>�   C2��Fݟ�Br�`�  C2�	<._Br���  C2�_Bt�Br }d�  C2�a[>�Br,<�  C3�j�f�Br7�h`  C3)W�UBrC��@  C3>��[��BrO{l   C3R\qE�Br[:�   C3g��{"�Brf�o�  C3{P�VBrr���  C3�U�4�Br~ys�  C3�v�:�	Br�8��  C3�Ǳfr�Br��w`  C3��E�Br���@  C3�u���Br�w{   C3��2a�;Br�6�   C3����Br��~�  C4�`A5Brж �  C4!È>bBr�u��  C4/�ŷwLBr�5�  C3���Br��`  C2�H��|Br��@  C2�`��eBss�   C1X�G&=�Bs3   C0��b��Bs"��  C0l��
ǬBs.��  C0��iDBs:q��  C/Ƭ��BsF1�  C/w�i13BsQ�`  C/B�o�Bs]�@  C/�7���Bsio�   C/r2�a�Bsu/   C.�����Bs���  C.��+�Bs���  C.�x�˅Bs�m��  C.�J&���Bs�-"�  C/��)Bs��`  C/|�o3Bs��&@  C/J���[eBs�k�   C/��ЈNBs�+*   C/Έ;D��Bs���  C0 ��x*�Bs�-�  C0]�m�Bs�i��  C05�OZ|�Bt)1�  C0O��#��Bt�`  C0iDw"ЪBt�5@  C0�H�?�}Bt%g�   C0��y�?Bt1'9   C0�6�s�HBt<��  C0����BtH�<�  C0�p[BtTe��  C1A3�'Bt`%@�  C1�RJ��Btk��`  C15o�<�~Btw�D@  C1Nu��bBt�c�   C1fƭ�vBt�#H   C1>UQi�Bt����  C1�i�LFKBt��K�  C1�g���Bt�a͠  C1�	PK2�Bt�!O�  C1�_����Bt���`  C1��>��BtՠS@  C2')��Bt�_�   C2!x8���Bt�W   C27�}��Bt����  C2N+�O%�Bu�Z�  C2dKo�Bu]ܠ  C2z;�a��Bu^�  C2�
�ٶ|Bu'��`  C2��G��Bu3�b@  C2�o��a�Bu?[�   C2����l�BuKf   C2����~�BuV���  C2��Ҍ�Bub�i�  C3�)�×BunY�  C3$��Buzm�  C39�nhBu���`  C3MAN�.Bu��q@  C3`PSl�Bu�W�   C3y�����Bu�u   C3����ģBu����  C3�o5/�Bu��x�  C3�J+�d�Bu�U��  C3�#��>Bu�|�  C3�ɒ�jBu���`  C3�PG���Bu@  C3�wʜ�Bu�T   C3�oַq�Bv�   C2�8'N�$Bv��  C2\d��Bv���  C1^>� �Bv*R	�  C0ծ��@Bv6��  C0l��2�BvA�`  C0!!S*B4BvM��@  C/�,���?BvYP   C/~���|�Bve�   C/?!���Bvp��  C/ �
Bv|���  C.�p�n�eBv�N�  C.��j׬�Bv���  C.�
�9�Bv��`  C.��� ��Bv���@  C.���$�OBv�L    C.���AzBv��   C.�5E�Bv��#�  C.��V��Bvڊ��  C.����Bv�J'�  C.�)�BGBv�	��  C/!�t��Bv��+`  C/L:ԫ?uBw	��@  C/~�H5��BwH/   C/���U�XBw!�   C/틚��aBw,�2�  C0���,Bw8���  C0,:�;�BwDF6�  C0F��E%
BwP��  C0`�\�cBw[�:`  C0z���g�Bwg��@  C0��� �BwsD>   C0���?Bw�   C0�LN\qBw��A�  C0�q, tBw����  C0��{8l6Bw�BE�  C1�d�"Bw�ǀ  C1,�[�SBw��I`  C1El5CdBwŀ�@  C1^�û��Bw�@M   C1w'N��)Bw���   C1�I3i�3Bw�P�  C1����LBw�~��  C1��K�%Bx >T�  C1�d=޵nBx�ր  C1�+=#��Bx�X`  C1��W�6Bx#|�@  C2�Z(�*Bx/<\   C2-�{�Bx:��   C2DQ��rBxF�_�  C2Z�6BxRz��  C2o4�ΪSBx^:c�  C2�yR�n�Bxi��  C2��)�NBxu�g`  C2�<y�UBx�x�@  C2�q�*PBx�8k   C2�Zx��Bx���   C2��w-�Bx��n�  C3�۪RaBx�v��  C39vtBx�6r�  C3+b.��Bx���  C3?��l��Bxӵv`  C3Sz���Bx�t�@  C3g����ZBx�4z   C3{���Bx���   C3�w��kBy�}�  C3�j:�7�Byr��  C3�f��(By2��  C3��'s��By%��  C3�7�CT�By1��`  C3�&(vCBy=q@  C3���+�ByI0�   C4֨�(ZByT�   C4 ��c�By`���  C48��ጹBylo�  C4P�B���Byx.��  C4\�zW�By���  C4oU�By���`  C4�ɏ�p�By�m@  C4��2���By�,�   C4���#By��   C4� ���qBy����  C3k�!w��By�k�  C2W��By�*��  C1xnUaL+By��!�  C0۠�'F]By���`  C0oߏ��By�i%@  C0�L:V�Bz(�   C/����Bz�)   C/z�	ZBz���  C/C���mBz(g,�  C/`�4��Bz4&��  C/E&�;Bz?�0�  C.��c�BzK��`  C.��<18VBzWe4@  C.���VcBzc$�   C.�j�	��Bzn�8   C/{�^��Bzz���  C/(wXh�^Bz�c;�  C/j��%��Bz�"��  C/�_s_Bz��?�  C/汥���