CDF  ?   
      time       lon       lat          	   CDI       @Climate Data Interface version 2.0.3 (https://mpimet.mpg.de/cdi)   Conventions       CF-1.5     source        #PISM development no-version-control    command      / /home/m/m300792/ModelCodes/mPISM_build/install_03-10-22/bin/pismr -i pism_-065000/pism_-065000.nc -stress_balance ssa+sia -y 100 -ys -65000 -o pism_-064900/pism_-064900.nc -surface_given_file forcing/smb_-064900.nc -surface given -ocean_th_file forcing/oce_-064900.nc -ocean th -topg_file pism_-065000/pism_-065000_new_topg.nc -bed_def prescribed -uplift_file pism_-065000/pism_-065000_dbdt.nc -extra_times 10 -extra_vars thk,usurf,bheatflx,velsurf_mag,velbase_mag,tillwat,bmelt,velbase,flux_divergence,velsurf,velbar,wvelsurf,wvelbase,climatic_mass_balance_cumulative,potential_climatic_mass_balance_cumulative,tempbase,tempsurf,mask,temppabase,tempicethk_basal,topg,nonneg_flux_cumulative,grounded_basal_flux_cumulative,floating_basal_flux_cumulative,discharge_flux_cumulative,Href,taud_mag,lat,lon,calving_mask,sea_level -extra_file pism_-064900/pism_-064900_extra -extra_split -ts_times yearly -ts_file pism_-064900/pism_-064900_ts.nc -topg_to_phi 15.0,30.0,-300.0,400.0 -ssafd_ksp_rtol 1e-7 -enthalpy_temperate_conductivity_ratio 0.001 -pseudo_plastic -pseudo_plastic_q .25 -pseudo_plastic_uthreshold 70 -till_effective_fraction_overburden 0.04 -tauc_slippery_grounding_lines -sia_e 8 -ssa_e 1 -ssa_upwind_factor 0.7 -calving eigen_calving,thickness_calving -thickness_calving_threshold 100 -part_redist -part_grid -cfbc -kill_icebergs -eigen_calving_K 1e16 -verbose 2 -sliding_scale_factor_reduces_tauc 1 -th_heat_coefficient 4e-5 -th_salinity_coefficient 4e-8 -o_format netcdf3 -sediment_file till_mask.nc -no_divide_tillphi -periodicity none -ocean_th_period 1 -o_order yxz
    proj4         j+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs     history_of_appended_files        �Fri Jul  8 17:11:28 2022: Appended file outFinal2 had following "history" attribute:
Fri Jul  8 17:11:27 2022: ncap2 -O -s lon_bnds=float(lon_bnds) outFinal1 outFinal2
Fri Jul  8 17:11:26 2022: ncap2 -O -s lat_bnds=float(lat_bnds) outFinal outFinal1
Fri Jul 08 17:11:25 2022: cdo -add thkNew thkOld outFinal
Fri Jul 08 17:11:24 2022: cdo -mul -selvar,thk pism_-064900/pism_-064900.nc MaskLaurentide thkNew
     NCO       _netCDF Operators version 5.0.6 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)    history       |Thu Sep 15 07:16:55 2022: cdo -mergetime File000001.nc File000002.nc File000003.nc File000004.nc File000005.nc File000006.nc File000007.nc File000008.nc File000009.nc File000010.nc File000011.nc File000012.nc File000013.nc File000014.nc File000015.nc File000016.nc File000017.nc File000018.nc File000019.nc File000020.nc File000021.nc File000022.nc File000023.nc File000024.nc File000025.nc File000026.nc File000027.nc File000028.nc File000029.nc File000030.nc File000031.nc File000032.nc File000033.nc File000034.nc File000035.nc File000036.nc File000037.nc File000038.nc File000039.nc File000040.nc File000041.nc File000042.nc File000043.nc File000044.nc File000045.nc File000046.nc File000047.nc File000048.nc File000049.nc File000050.nc File000051.nc File000052.nc File000053.nc File000054.nc File000055.nc File000056.nc File000057.nc File000058.nc File000059.nc File000060.nc File000061.nc File000062.nc File000063.nc File000064.nc File000065.nc File000066.nc File000067.nc File000068.nc File000069.nc File000070.nc File000071.nc File000072.nc File000073.nc File000074.nc File000075.nc File000076.nc File000077.nc File000078.nc File000079.nc File000080.nc File000081.nc File000082.nc File000083.nc File000084.nc File000085.nc File000086.nc File000087.nc File000088.nc File000089.nc File000090.nc File000091.nc File000092.nc File000093.nc File000094.nc File000095.nc File000096.nc File000097.nc File000098.nc File000099.nc File000100.nc File000101.nc File000102.nc File000103.nc File000104.nc File000105.nc File000106.nc File000107.nc File000108.nc File000109.nc File000110.nc File000111.nc File000112.nc File000113.nc File000114.nc File000115.nc File000116.nc File000117.nc File000118.nc File000119.nc File000120.nc File000121.nc File000122.nc File000123.nc File000124.nc File000125.nc File000126.nc File000127.nc File000128.nc File000129.nc File000130.nc File000131.nc File000132.nc File000133.nc File000134.nc File000135.nc File000136.nc File000137.nc File000138.nc File000139.nc File000140.nc File000141.nc File000142.nc File000143.nc File000144.nc File000145.nc File000146.nc File000147.nc File000148.nc File000149.nc File000150.nc File000151.nc File000152.nc File000153.nc File000154.nc File000155.nc File000156.nc File000157.nc File000158.nc File000159.nc File000160.nc File000161.nc File000162.nc File000163.nc File000164.nc File000165.nc File000166.nc File000167.nc File000168.nc File000169.nc File000170.nc File000171.nc File000172.nc File000173.nc File000174.nc File000175.nc File000176.nc File000177.nc File000178.nc File000179.nc File000180.nc File000181.nc File000182.nc File000183.nc File000184.nc File000185.nc File000186.nc File000187.nc File000188.nc File000189.nc File000190.nc File000191.nc File000192.nc File000193.nc File000194.nc File000195.nc File000196.nc File000197.nc File000198.nc File000199.nc File000200.nc File000201.nc File000202.nc File000203.nc File000204.nc File000205.nc File000206.nc File000207.nc File000208.nc File000209.nc File000210.nc File000211.nc File000212.nc File000213.nc File000214.nc File000215.nc File000216.nc File000217.nc File000218.nc File000219.nc File000220.nc File000221.nc File000222.nc File000223.nc File000224.nc File000225.nc File000226.nc File000227.nc File000228.nc File000229.nc File000230.nc File000231.nc File000232.nc File000233.nc File000234.nc File000235.nc File000236.nc File000237.nc File000238.nc File000239.nc File000240.nc File000241.nc File000242.nc File000243.nc File000244.nc File000245.nc File000246.nc File000247.nc File000248.nc File000249.nc File000250.nc File000251.nc File000252.nc File000253.nc File000254.nc File000255.nc File000256.nc File000257.nc File000258.nc File000259.nc File000260.nc File000261.nc File000262.nc File000263.nc File000264.nc File000265.nc File000266.nc File000267.nc File000268.nc File000269.nc File000270.nc File000271.nc File000272.nc File000273.nc File000274.nc File000275.nc File000276.nc File000277.nc File000278.nc File000279.nc File000280.nc File000281.nc File000282.nc File000283.nc File000284.nc File000285.nc File000286.nc File000287.nc File000288.nc File000289.nc File000290.nc File000291.nc File000292.nc File000293.nc File000294.nc File000295.nc File000296.nc File000297.nc File000298.nc File000299.nc File000300.nc File000301.nc File000302.nc File000303.nc File000304.nc File000305.nc File000306.nc File000307.nc File000308.nc File000309.nc File000310.nc File000311.nc File000312.nc File000313.nc File000314.nc File000315.nc File000316.nc File000317.nc File000318.nc File000319.nc File000320.nc File000321.nc File000322.nc File000323.nc File000324.nc File000325.nc File000326.nc File000327.nc File000328.nc File000329.nc File000330.nc File000331.nc File000332.nc File000333.nc File000334.nc File000335.nc File000336.nc File000337.nc File000338.nc File000339.nc File000340.nc File000341.nc File000342.nc File000343.nc File000344.nc File000345.nc File000346.nc File000347.nc File000348.nc File000349.nc File000350.nc File000351.nc File000352.nc File000353.nc File000354.nc File000355.nc File000356.nc File000357.nc File000358.nc File000359.nc File000360.nc File000361.nc File000362.nc File000363.nc File000364.nc File000365.nc File000366.nc File000367.nc File000368.nc File000369.nc File000370.nc File000371.nc File000372.nc File000373.nc File000374.nc File000375.nc File000376.nc File000377.nc File000378.nc File000379.nc File000380.nc File000381.nc File000382.nc File000383.nc File000384.nc File000385.nc File000386.nc File000387.nc File000388.nc File000389.nc File000390.nc File000391.nc File000392.nc File000393.nc File000394.nc File000395.nc File000396.nc File000397.nc File000398.nc File000399.nc File000400.nc File000401.nc File000402.nc File000403.nc File000404.nc File000405.nc File000406.nc File000407.nc File000408.nc File000409.nc File000410.nc File000411.nc File000412.nc File000413.nc File000414.nc File000415.nc File000416.nc File000417.nc File000418.nc File000419.nc File000420.nc File000421.nc File000422.nc File000423.nc File000424.nc File000425.nc File000426.nc File000427.nc File000428.nc File000429.nc File000430.nc File000431.nc File000432.nc File000433.nc File000434.nc File000435.nc File000436.nc File000437.nc File000438.nc File000439.nc File000440.nc File000441.nc File000442.nc File000443.nc File000444.nc File000445.nc File000446.nc File000447.nc File000448.nc File000449.nc File000450.nc File000451.nc File000452.nc File000453.nc File000454.nc File000455.nc File000456.nc File000457.nc File000458.nc File000459.nc File000460.nc File000461.nc File000462.nc File000463.nc File000464.nc File000465.nc File000466.nc File000467.nc File000468.nc File000469.nc File000470.nc File000471.nc File000472.nc File000473.nc File000474.nc File000475.nc File000476.nc File000477.nc File000478.nc File000479.nc File000480.nc File000481.nc File000482.nc File000483.nc File000484.nc File000485.nc File000486.nc File000487.nc File000488.nc File000489.nc File000490.nc File000491.nc File000492.nc File000493.nc File000494.nc File000495.nc File000496.nc File000497.nc File000498.nc File000499.nc File000500.nc File000501.nc File000502.nc File000503.nc File000504.nc File000505.nc File000506.nc File000507.nc File000508.nc File000509.nc File000510.nc File000511.nc File000512.nc File000513.nc File000514.nc File000515.nc File000516.nc File000517.nc File000518.nc File000519.nc File000520.nc File000521.nc File000522.nc File000523.nc File000524.nc File000525.nc File000526.nc File000527.nc File000528.nc File000529.nc File000530.nc File000531.nc File000532.nc File000533.nc File000534.nc File000535.nc File000536.nc File000537.nc File000538.nc File000539.nc File000540.nc File000541.nc File000542.nc File000543.nc File000544.nc File000545.nc File000546.nc File000547.nc File000548.nc File000549.nc File000550.nc File000551.nc File000552.nc File000553.nc File000554.nc File000555.nc File000556.nc File000557.nc File000558.nc File000559.nc File000560.nc File000561.nc File000562.nc File000563.nc File000564.nc File000565.nc File000566.nc File000567.nc File000568.nc File000569.nc File000570.nc File000571.nc File000572.nc File000573.nc File000574.nc File000575.nc IceVolumeHudson.nc
Thu Sep 15 07:09:26 2022: cdo -fldsum -mulc,10000 -mulc,10000 -selindexbox,279,477,230,362 -selvar,thk /work/ba0989/m300792/HE_Runs/ExperimentsComposite//HE77//pism_-064900/pism_-064900.nc Tmp/File000001.nc   CDO       @Climate Data Operators version 2.0.3 (https://mpimet.mpg.de/cdo)         time                standard_name         time   	long_name         time   units         seconds since 1-1-1    calendar      365_day    axis      T               -�   lon                standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X               -�   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y               -�   thk                       standard_name         land_ice_thickness     	long_name         land ice thickness     units         m      pism_intent       model_state             -�                �}Ȁ@�  C4�(w-o�}���   C4���v��}�=   C4�ÚB}��}�A�@  C3��w��L�}��9`  C2�4��}�·�  C1�x �p�}�5�  C1A��kak�}vC��  C0�w�A�}j�1�  C0���k�}^İ   C0��s*���}S.   C0��D��}GE�@  C0}�(�g�};�*`  C0�D]I���}/ƨ�  C0�)�*b�}$&�  C0���܀g�}G��  C0���#jn�}�"�  C0ǈ(��>�} ȡ   C0�~��Ew�|�	   C0�R de�|�I�@  C0����v��|݊`  C0�=���b�|�ʙ�  C1�0�f��|��  C1E]^Y�|�K��  C1p��|���  C1�.~)���|�̒   C1�,�}P��|�   C1�Ǩp�|�M�@  C2�2B��|�`  C2*3�J!^�|sΊ�  C2D�^!
��|h�  C2]?Nb��|\O��  C2p4�,)�|P��  C2�)9��|DЃ   C2�����|9   C2�۴Bo��|-Q@  C2��Yj���|!��`  C2�����|�{�  C2�蘼cn�|
��  C2�3�:p�{�Sw�  C2�lf"�(�{���  C2����\��{��t   C2�O	�ޑ�{��   C3w_D,��{�Up@  C327x|��{Õ�`  C3Y���L�{��l�  C3|��lO��{��  C3��Y��^�{�Wh�  C3������{����  C3�=)Ɲ*�{��e   C3厈*v�{}�   C3�֥v�{qYa@  C4	3(XK��{e��`  C4C�]�{Y�]�  C4-UrO�{N۠  C4��W��{B[Y�  C4B^�u�{6���  C4�.�ѥ�{*�V   C4�G�{�   C4�g_�:�{]R@  C4ð!��{��`  C4*̭�x��z��N�  C4M%�e��z�̠  C4s�K�z�_J�  C4� �X	��z؟��  C4�nFV(U�z��G   C3�V���z� �   C2�#B=�z�aC@  C1h.���k�z���`  C0���o?z�z��?�  C0���e��z�"��  C0���WG�z�c;�  C0iˊ�Z�zz���  C0T`��i��zn�8   C0F�[�x(�zc$�   C0U8���zWe4@  C0hAJ�:m�zK��`  C0uI�8�H�z?�0�  C0~]�s�z4&��  C0�N>���z(g,�  C0�T�M[��z���  C0������z�)   C0���e��z(�   C0������y�i%@  C1'������y���`  C1S�'nF�y��!�  C1}Sgf���y�*��  C1��ԭ���y�k�  C1Ø�Q��y����  C1������y��   C28XV���y�,�   C2I�{`��y�m@  C23+�;H	�y���`  C2<����=�y���  C2G9�|�yx.��  C2O�����ylo�  C2T��*���y`���  C2W��o`z�yT�   C2X1���@�yI0�   C2YR�|�y=q@  C2^�^h���y1��`  C2t(@��n�y%��  C2��� ��y2��  C2Ɔl�)�yr��  C2�^�gg��y�}�  C3�S�Z^�x���   C3?Cy{��x�4z   C3`���ad�x�t�@  C3~�3��[�xӵv`  C3��̝+�x���  C3��f2`�x�6r�  C3�O�m;�x�v��  C3�ɵ@�2�x��n�  C3�$;���x���   C3���T{�x�8k   C3������x�x�@  C3�kP:��xu�g`  C3�N�@�xi��  C3��	 ��x^:c�  C3�?�BA�xRz��  C3�yi=��xF�_�  C3�e0�]C�x:��   C4!���P.�x/<\   C4Ha
AJ�x#|�@  C4l�3c~�x�X`  C4���W�x�ր  C3+���9�x >T�  C21Kh�G��w�~��  C1gIV�w�P�  C0�J�4�w���   C0����y �w�@M   C0i�X���wŀ�@  C0F�Q�g��w��I`  C01�~�0�w�ǀ  C0�[
j��w�BE�  C0;�:�w����  C01����w��A�  C0@���e�w�   C0Kb��F�wsD>   C0Ty"\��wg��@  C0]�̫v��w[�:`  C0s�e�� �wP��  C0��"��wDF6�  C0̕Ŕ��w8���  C0��U����w,�2�  C1'��k9�w!�   C1Qq�
O~�wH/   C1x���j�w	��@  C1�"l��v��+`  C1��"W�,�v�	��  C1֘��,��v�J'�  C1�RH���vڊ��  C2
�{�� �v��#�  C2	܆vT�v��   C2!�6<3)�v�L    C2+	����v���@  C21,P�3�v��`  C24l;�a?�v���  C26_%����v�N�  C2:� Չ�v|���  C2E�h���vp��  C2_ �,#Q�ve�   C2{�jN���vYP   C2�_9f�vM��@  C2�֨���vA�`  C2�
��v�v6��  C3j`9ML�v*R	�  C3@�5o�v���  C3\���n�v��  C3v���$��v�   C3�/�j(�u�T   C3�P�ϛ~�u@  C3��P���u���`  C3�T}���u�|�  C3�)�{u��u�U��  C3ȅ�@���u��x�  C3���Q�e�u����  C3��Uī��u�u   C3�A%��A�u�W�   C3���H�u��q@  C3�#*d/;�u���`  C3׾8g��uzm�  C3��\�k��unY�  C4!�X�7��ub�i�  C4F��m*�uV���  C4g5d�uKf   C3m C���u?[�   C2^��T�u3�b@  C1|�fBMF�u'��`  C0�	�e���u^�  C0��C�T�u]ܠ  C0E_�����u�Z�  C0m�V�-�t����  C/�8���B�t�W   C/�
�C+_�t�_�   C/���F���tՠS@  C/z�lԣ��t���`  C/y��i��t�!O�  C/��k��q�t�a͠  C/��l��w�t��K�  C/�X���t����  C0����t�#H   C0/�W���t�c�   C0[�����tw�D@  C0���cU��tk��`  C0�23��t`%@�  C0���9�J�tTe��  C1	�%��tH�<�  C1,����t<��  C1J����^�t1'9   C1e��_��t%g�   C1~1�u���t�5@  C1�����t�`  C1���s��t)1�  C1��hO���s�i��  C1���g�h�s�-�  C1�<�C`B�s���  C1��]��6�s�+*   C1�ٮ>��s�k�   C1ǣÇT��s��&@  C1�i�ڵ��s��`  C1�4r��P�s�-"�  C2����s�m��  C29H!Qd��s���  C2b�q�:&�s���  C2�.UI�2�su/   C2�Q����sio�   C2� ��S�s]�@  C2�ܟp�]�sQ�`  C3G��V��sF1�  C3+WY�4�s:q��  C3B^$��s.��  C3W\9�_�s"��  C3`Y
����s3   C3i���ss�   C3p
�����r��@  C3sVn�m�r��`  C3t�.W��r�5�  C3s*�'@��r�u��  C3r���rж �  C3u�K���r��~�  C3�A}bSl�r�6�   C3�h`��l�r�w{   C3�NE�5�r���@  C3��j5>�r��w`  C4��_P��r�8��  C46��YE��r~ys�  C3چ��
�rr���  C2.k����rf�o�  C1V��Mrf�r[:�   C0ŋ��j�rO{l   C0d�+-��rC��@  C0&���K�r7�h`  C/�6����r,<�  C/��HM��r }d�  C/�W��L�r���  C/`�| 0��r�`�  C/JE�p��q�>�   C/j��}��q�]   C/�6KDv��q��@  C/�:�3���q� Y`  C/�>����q�@׀  C09ĩ��qU�  C0K�77~^�q����  C0{Pa��q�Q�  C0���t��q�B�   C0��r{��q��N   C0���#��q���@  C1��)�q|J`  C1@�S��i�qpDȀ  C1\�K����qd�F�  C1v�Υfq�qX���  C1��,��qMB�  C1�1_JlM�qAF�   C1�����q5�?   C1���0�q)ǽ@  C1��_��q;`  C1�4�I'��qH��  C1�@ݽ!�q�7�  C1�,�Y���p�ɵ�  C1аXw��p�
3�  C1䮷�#��p�J�   C2@r'`�p׋0   C26~�܋0�p�ˮ@  C2a����p�,`  C2���͠�p�L��  C2�yO���p��(�  C2�<+#E�p�ͦ�  C2��~��p�$�  C3���h��p�N�   C3&�f����py�!   C3:���j�pmϟ@  C3G�l>��pb`  C3Uh"�a��pVP��  C3]����Z�pJ��  C3coI_kM�p>ї�  C3f�Z!J�p3�  C3g_�K=��p'R�   C3f,�H:�p�   C3e5���U�pӐ@  C3i)8�I�p`  C3{�Dֈ8�o�   C3���X���o�*@  C3���͐�o���  C3�\|��o�,�  C4Ex���o��
   C3Ձb�o{.@  C2�F0��oc��  C1��/ƣ��oL/��  C1��"��o4��   C0��a�/�o1�@  C04@�(k��o��  C/�ʏ��S�n�3��  C/�cތ���nִ�   C/rX���n�5�@  C/E oP�n���  C/$��~�+�n�7��  C/�&$�&�nx��   C/ 8;9��na9�@  C/</�BI�nI�Հ  C/]���%X�n2;��  C/�\�2ڲ�n��   C/�u!(���n=�@  C0"+��2�m�ƀ  C0Q��b�m�?��  C0+5p���m���   C0��Z���m�A�@  C0������m�·�  C0�j��x��mvC��  C1fS��m^İ   C11����mGE�@  C1HfV���m/ƨ�  C1]�H��N�mG��  C1q��V�m ȡ   C1�5�EP,�l�I�@  C1����'�l�ʙ�  C1��*�'��l�K��  C1�e����l�̒   C1�j_<m�l�M�@  C1������lsΊ�  C1��v|���l\O��  C1�'�h���lDЃ   C1���Ȇ��l-Q@  C2	����l�{�  C24QBM�@�k�Sw�  C2[d6� ��k��t   C2�ጉW��k�Up@  C2�-IQ]z�k��l�  C2��y����k�Wh�  C2�{�d���k��e   C2�7�Y�kqYa@  C3Iz���kY�]�  C3 
3ḍ�kB[Y�  C3.�~oΰ�k*�V   C398�Ԯ�k]R@  C3@�ģ�s�j��N�  C3E#��c�j�_J�  C3Fbkm�.�j��G   C3G��-#��j�aC@  C3J��f�B�j��?�  C3R��o�6�j�c;�  C3l0���jn�8   C3�����^�jWe4@  C3�y�?&�j?�0�  C3τK��
�j(g,�  C3�d�%���j�)   C4�/*E�i�i%@  C46���]M�i��!�  C4S��MJ��i�k�  C4kŏt^��i��   C4}�q��D�i�m@  C4+C5��i���  C3<�ݚX�ilo�  C2/�b�3�iT�   C1v�|�Z�i=q@  C0���a��i%��  C0���F��ir��  C0U�#�g��h���   C0&C$���h�t�@  C0�iq���h���  C/�{�+���h�v��  C/����7I�h���   C/rd7�w��h�x�@  C/ZBKK��hi��  C/K�c���hRz��  C/E_����h:��   C/@�|ش��h#|�@  C/3�v�D�h�ր  C/P;	���g�~��  C/LzԹ�g���   C.��w�*��gŀ�@  C.ۊWI�
�g�ǀ  C.��#��g����  C.�j�Q<��g�   C/'s���gg��@  C/E�����gP��  C/T$�����g8���  C/dLa���g!�   C/sa���4�g	��@  C/������f�	��  C/�]���s�fڊ��  C/����.��f��   C0
��ٿ��f���@  C0:�5��f���  C0j�%���f|���  C0�}m�V�fe�   C0��T
(��fM��@  C0�z}s֡�f6��  C1��3�f���  C1,
i�f�   C1G(~��e@  C1`M�����e�|�  C1q�����e��x�  C1����&��e�u   C1�WFI���e��q@  C1��=��u�ezm�  C1�_uvO��eb�i�  C1�>+զ��eKf   C1��+J�e3�b@  C1��؄ߦ�e^�  C1�����n�e�Z�  C1�Z���<�d�W   C1� Ύ���dՠS@  C2�-�d�!O�  C2FW�k�-�d��K�  C2o��ҝ�d�#H   C2��r]���dw�D@  C2��'U]��d`%@�  C2�+��k��dH�<�  C2������d1'9   C3:�ή��d�5@  C3$�s
�d)1�  C35�����c�-�  C3B��P��c�+*   C3LDc�c��&@  C3R`˺��c�-"�  C3X��/ �c���  C3]��4���cu/   C3_�?��:�c]�@  C3[׫��L�cF1�  C3[E� �c.��  C3m�+}�c3   C3��ǌ���b��@  C3��\���b�5�  C3،��E��bж �  C3�n#o���b�6�   C3BGނ�g�b���@  C2D��<��b�8��  C1�;B+B�br���  C0�TD���b[:�   C0b��nS�bC��@  C0����b,<�  C/�9|��b���  C/hy��@^�a�>�   C/&��޳��a��@  C.�~���a�@׀  C.�;J�S��a����  C.��G�a�B�   C.����8K�a���@  C/������apDȀ  C/*U.H,�aX���  C/m)�H� �aAF�   C/��ˇ��a)ǽ@  C06my�i�aH��  C0=p�n���`�ɵ�  C0k��<��`�J�   C0��Y��T�`�ˮ@  C0��O����`�L��  C0��ZG��`�ͦ�  C1��^�2�`�N�   C1!m$�]�`mϟ@  C1:��R�E�`VP��  C1Sy�Hn�`>ї�  C1i�@�k��`'R�   C1v%4U��`Ӑ@  C1|����r�_�   C1�!�����_���  C1���'��_��
   C1���1�_c��  C1��K�D��_4��   C1��7�g�_��  C1�D��R�^ִ�   C1֐�&.�^���  C2�d��Q�^x��   C2-3(�D��^I�Հ  C2Vc�ݬ��^��   C2|�rAh��]�ƀ  C2��뫑��]���   C2��Y���]�·�  C2�pK�A��]^İ   C2���[��]/ƨ�  C3	�:aq�] ȡ   C3�����\�ʙ�  C3-�!a���\�̒   C3:�e˹G�\sΊ�  C39J�o��\DЃ   C3;�����\�{�  C3:�1��[��t   C39R%����[��l�  C380K�;�[��e   C3;ӓi�J�[Y�]�  C3N�ZN���[*�V   C3r�	���Z��N�  C3���^V�Z��G   C3����L(�Z��?�  C3��g(��Zn�8   C3��]�u�Z?�0�  C2�C\&��Z�)   C1��Y��Y��!�  C0���\��Y��   C0o�'���Y���  C0Ǉ��V�YT�   C/����Y%��  C/i����U�X���   C/s�8���X���  C.���`f�X���   C.р�(R}�Xi��  C.�ovb�q�X:��   C.�Zn�L��X�ր  C.���VŎ�W���   C/B9!q��W�ǀ  C/G���2�W�   C/�eA���WP��  C/���H��W!�   C0%�G��5�V�	��  C0T�^�/�V��   C0��í��V���  C0�آf���Ve�   C0�,��\��V6��  C0��"�<J�V�   C1�Y��1�U�|�  C1�PWʷ�U�u   C10,���Uzm�  C1A������UKf   C1N��{r?�U^�  C1Yҫ��O�T�W   C1`�_���T�!O�  C1d��\��T�#H   C1j��z��T`%@�  C1p6�i8�T1'9   C1z��J��T)1�  C1��6��7�S�+*   C1��|����S�-"�  C1ں��u��Su/   C2�{y���SF1�  C20i�f��S3   C2V�M��f�R�5�  C2z��A�F�R�6�   C2��#i��R�8��  C2��&���R[:�   C2�
�G��R,<�  C2�i���R�Q�>�   C2��~;B�Q�@׀  C3~� R�Q�B�   C3H"�}�QpDȀ  C3]�Ȉ��QAF�   C3����x�QH��  C3 7H;���P�J�   C3�'NL�P�L��  C3PC�_d�P�N�   C3#9��#��PVP��  C36� �/��P'R�   C3Z�fJ/��O�   C3�Ol�ح�O��
   C3���ʼ�O4��   C3�����Nִ�   C2t����+�Nx��   C1��S��F�N��   C0��C7y��M���   C0[<��4>�M^İ   C0zH����M ȡ   C/�_�.[�L�̒   C/'���LDЃ   C.�ӊ�+��K��t   C.΁���o�K��e   C.�+��d