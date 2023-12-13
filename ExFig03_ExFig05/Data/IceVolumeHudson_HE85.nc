CDF  ?   
      time       lon       lat          	   CDI       @Climate Data Interface version 2.0.3 (https://mpimet.mpg.de/cdi)   Conventions       CF-1.5     source        #PISM development no-version-control    command      / /home/m/m300792/ModelCodes/mPISM_build/install_03-10-22/bin/pismr -i pism_-065000/pism_-065000.nc -stress_balance ssa+sia -y 100 -ys -65000 -o pism_-064900/pism_-064900.nc -surface_given_file forcing/smb_-064900.nc -surface given -ocean_th_file forcing/oce_-064900.nc -ocean th -topg_file pism_-065000/pism_-065000_new_topg.nc -bed_def prescribed -uplift_file pism_-065000/pism_-065000_dbdt.nc -extra_times 10 -extra_vars thk,usurf,bheatflx,velsurf_mag,velbase_mag,tillwat,bmelt,velbase,flux_divergence,velsurf,velbar,wvelsurf,wvelbase,climatic_mass_balance_cumulative,potential_climatic_mass_balance_cumulative,tempbase,tempsurf,mask,temppabase,tempicethk_basal,topg,nonneg_flux_cumulative,grounded_basal_flux_cumulative,floating_basal_flux_cumulative,discharge_flux_cumulative,Href,taud_mag,lat,lon,calving_mask,sea_level -extra_file pism_-064900/pism_-064900_extra -extra_split -ts_times yearly -ts_file pism_-064900/pism_-064900_ts.nc -topg_to_phi 15.0,30.0,-300.0,400.0 -ssafd_ksp_rtol 1e-7 -enthalpy_temperate_conductivity_ratio 0.001 -pseudo_plastic -pseudo_plastic_q .25 -pseudo_plastic_uthreshold 70 -till_effective_fraction_overburden 0.04 -tauc_slippery_grounding_lines -sia_e 8 -ssa_e 1 -ssa_upwind_factor 0.7 -calving eigen_calving,thickness_calving -thickness_calving_threshold 100 -part_redist -part_grid -cfbc -kill_icebergs -eigen_calving_K 1e16 -verbose 2 -sliding_scale_factor_reduces_tauc 1 -th_heat_coefficient 4e-5 -th_salinity_coefficient 4e-8 -o_format netcdf3 -sediment_file till_mask.nc -no_divide_tillphi -periodicity none -ocean_th_period 1 -o_order yxz
    proj4         j+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs     history_of_appended_files        �Wed Jul 20 09:54:51 2022: Appended file outFinal2 had following "history" attribute:
Wed Jul 20 09:54:50 2022: ncap2 -O -s lon_bnds=float(lon_bnds) outFinal1 outFinal2
Wed Jul 20 09:54:49 2022: ncap2 -O -s lat_bnds=float(lat_bnds) outFinal outFinal1
Wed Jul 20 09:54:48 2022: cdo -add thkNew thkOld outFinal
Wed Jul 20 09:54:47 2022: cdo -mul -selvar,thk pism_-064900/pism_-064900.nc MaskLaurentide thkNew
     NCO       _netCDF Operators version 5.0.6 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)    history       |Thu Sep 15 07:21:50 2022: cdo -mergetime File000001.nc File000002.nc File000003.nc File000004.nc File000005.nc File000006.nc File000007.nc File000008.nc File000009.nc File000010.nc File000011.nc File000012.nc File000013.nc File000014.nc File000015.nc File000016.nc File000017.nc File000018.nc File000019.nc File000020.nc File000021.nc File000022.nc File000023.nc File000024.nc File000025.nc File000026.nc File000027.nc File000028.nc File000029.nc File000030.nc File000031.nc File000032.nc File000033.nc File000034.nc File000035.nc File000036.nc File000037.nc File000038.nc File000039.nc File000040.nc File000041.nc File000042.nc File000043.nc File000044.nc File000045.nc File000046.nc File000047.nc File000048.nc File000049.nc File000050.nc File000051.nc File000052.nc File000053.nc File000054.nc File000055.nc File000056.nc File000057.nc File000058.nc File000059.nc File000060.nc File000061.nc File000062.nc File000063.nc File000064.nc File000065.nc File000066.nc File000067.nc File000068.nc File000069.nc File000070.nc File000071.nc File000072.nc File000073.nc File000074.nc File000075.nc File000076.nc File000077.nc File000078.nc File000079.nc File000080.nc File000081.nc File000082.nc File000083.nc File000084.nc File000085.nc File000086.nc File000087.nc File000088.nc File000089.nc File000090.nc File000091.nc File000092.nc File000093.nc File000094.nc File000095.nc File000096.nc File000097.nc File000098.nc File000099.nc File000100.nc File000101.nc File000102.nc File000103.nc File000104.nc File000105.nc File000106.nc File000107.nc File000108.nc File000109.nc File000110.nc File000111.nc File000112.nc File000113.nc File000114.nc File000115.nc File000116.nc File000117.nc File000118.nc File000119.nc File000120.nc File000121.nc File000122.nc File000123.nc File000124.nc File000125.nc File000126.nc File000127.nc File000128.nc File000129.nc File000130.nc File000131.nc File000132.nc File000133.nc File000134.nc File000135.nc File000136.nc File000137.nc File000138.nc File000139.nc File000140.nc File000141.nc File000142.nc File000143.nc File000144.nc File000145.nc File000146.nc File000147.nc File000148.nc File000149.nc File000150.nc File000151.nc File000152.nc File000153.nc File000154.nc File000155.nc File000156.nc File000157.nc File000158.nc File000159.nc File000160.nc File000161.nc File000162.nc File000163.nc File000164.nc File000165.nc File000166.nc File000167.nc File000168.nc File000169.nc File000170.nc File000171.nc File000172.nc File000173.nc File000174.nc File000175.nc File000176.nc File000177.nc File000178.nc File000179.nc File000180.nc File000181.nc File000182.nc File000183.nc File000184.nc File000185.nc File000186.nc File000187.nc File000188.nc File000189.nc File000190.nc File000191.nc File000192.nc File000193.nc File000194.nc File000195.nc File000196.nc File000197.nc File000198.nc File000199.nc File000200.nc File000201.nc File000202.nc File000203.nc File000204.nc File000205.nc File000206.nc File000207.nc File000208.nc File000209.nc File000210.nc File000211.nc File000212.nc File000213.nc File000214.nc File000215.nc File000216.nc File000217.nc File000218.nc File000219.nc File000220.nc File000221.nc File000222.nc File000223.nc File000224.nc File000225.nc File000226.nc File000227.nc File000228.nc File000229.nc File000230.nc File000231.nc File000232.nc File000233.nc File000234.nc File000235.nc File000236.nc File000237.nc File000238.nc File000239.nc File000240.nc File000241.nc File000242.nc File000243.nc File000244.nc File000245.nc File000246.nc File000247.nc File000248.nc File000249.nc File000250.nc File000251.nc File000252.nc File000253.nc File000254.nc File000255.nc File000256.nc File000257.nc File000258.nc File000259.nc File000260.nc File000261.nc File000262.nc File000263.nc File000264.nc File000265.nc File000266.nc File000267.nc File000268.nc File000269.nc File000270.nc File000271.nc File000272.nc File000273.nc File000274.nc File000275.nc File000276.nc File000277.nc File000278.nc File000279.nc File000280.nc File000281.nc File000282.nc File000283.nc File000284.nc File000285.nc File000286.nc File000287.nc File000288.nc File000289.nc File000290.nc File000291.nc File000292.nc File000293.nc File000294.nc File000295.nc File000296.nc File000297.nc File000298.nc File000299.nc File000300.nc File000301.nc File000302.nc File000303.nc File000304.nc File000305.nc File000306.nc File000307.nc File000308.nc File000309.nc File000310.nc File000311.nc File000312.nc File000313.nc File000314.nc File000315.nc File000316.nc File000317.nc File000318.nc File000319.nc File000320.nc File000321.nc File000322.nc File000323.nc File000324.nc File000325.nc File000326.nc File000327.nc File000328.nc File000329.nc File000330.nc File000331.nc File000332.nc File000333.nc File000334.nc File000335.nc File000336.nc File000337.nc File000338.nc File000339.nc File000340.nc File000341.nc File000342.nc File000343.nc File000344.nc File000345.nc File000346.nc File000347.nc File000348.nc File000349.nc File000350.nc File000351.nc File000352.nc File000353.nc File000354.nc File000355.nc File000356.nc File000357.nc File000358.nc File000359.nc File000360.nc File000361.nc File000362.nc File000363.nc File000364.nc File000365.nc File000366.nc File000367.nc File000368.nc File000369.nc File000370.nc File000371.nc File000372.nc File000373.nc File000374.nc File000375.nc File000376.nc File000377.nc File000378.nc File000379.nc File000380.nc File000381.nc File000382.nc File000383.nc File000384.nc File000385.nc File000386.nc File000387.nc File000388.nc File000389.nc File000390.nc File000391.nc File000392.nc File000393.nc File000394.nc File000395.nc File000396.nc File000397.nc File000398.nc File000399.nc File000400.nc File000401.nc File000402.nc File000403.nc File000404.nc File000405.nc File000406.nc File000407.nc File000408.nc File000409.nc File000410.nc File000411.nc File000412.nc File000413.nc File000414.nc File000415.nc File000416.nc File000417.nc File000418.nc File000419.nc File000420.nc File000421.nc File000422.nc File000423.nc File000424.nc File000425.nc File000426.nc File000427.nc File000428.nc File000429.nc File000430.nc File000431.nc File000432.nc File000433.nc File000434.nc File000435.nc File000436.nc File000437.nc File000438.nc File000439.nc File000440.nc File000441.nc File000442.nc File000443.nc File000444.nc File000445.nc File000446.nc File000447.nc File000448.nc File000449.nc File000450.nc File000451.nc File000452.nc File000453.nc File000454.nc File000455.nc File000456.nc File000457.nc File000458.nc File000459.nc File000460.nc File000461.nc File000462.nc File000463.nc File000464.nc File000465.nc File000466.nc File000467.nc File000468.nc File000469.nc File000470.nc File000471.nc File000472.nc File000473.nc File000474.nc File000475.nc File000476.nc File000477.nc File000478.nc File000479.nc File000480.nc File000481.nc File000482.nc File000483.nc File000484.nc File000485.nc File000486.nc File000487.nc File000488.nc File000489.nc File000490.nc File000491.nc File000492.nc File000493.nc File000494.nc File000495.nc File000496.nc File000497.nc File000498.nc File000499.nc File000500.nc File000501.nc File000502.nc File000503.nc File000504.nc File000505.nc File000506.nc File000507.nc File000508.nc File000509.nc File000510.nc File000511.nc File000512.nc File000513.nc File000514.nc File000515.nc File000516.nc File000517.nc File000518.nc File000519.nc File000520.nc File000521.nc File000522.nc File000523.nc File000524.nc File000525.nc File000526.nc File000527.nc File000528.nc File000529.nc File000530.nc File000531.nc File000532.nc File000533.nc File000534.nc File000535.nc File000536.nc File000537.nc File000538.nc File000539.nc File000540.nc File000541.nc File000542.nc File000543.nc File000544.nc File000545.nc File000546.nc File000547.nc File000548.nc File000549.nc File000550.nc File000551.nc File000552.nc File000553.nc File000554.nc File000555.nc File000556.nc File000557.nc File000558.nc File000559.nc File000560.nc File000561.nc File000562.nc File000563.nc File000564.nc File000565.nc File000566.nc File000567.nc File000568.nc File000569.nc File000570.nc File000571.nc File000572.nc File000573.nc File000574.nc File000575.nc IceVolumeHudson.nc
Thu Sep 15 07:14:58 2022: cdo -fldsum -mulc,10000 -mulc,10000 -selindexbox,279,477,230,362 -selvar,thk /work/ba0989/m300792/HE_Runs/ExperimentsComposite//HE85//pism_-064900/pism_-064900.nc Tmp/File000001.nc   CDO       @Climate Data Operators version 2.0.3 (https://mpimet.mpg.de/cdo)         time                standard_name         time   	long_name         time   units         seconds since 1-1-1    calendar      365_day    axis      T               -�   lon                standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X               -�   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y               -�   thk                       standard_name         land_ice_thickness     	long_name         land ice thickness     units         m      pism_intent       model_state             -�                �}Ȁ@�  C4�(w-o�}���   C4� D���}�=   C4��b�S�}�A�@  C4�]t�W�}��9`  C4B�;���}�·�  C3 ����}�5�  C2waoπ�}vC��  C1���9�}j�1�  C1Z�Գ�}^İ   C1�ꟍ�}S.   C0��y��C�}GE�@  C0���*��};�*`  C0�i�l��}/ƨ�  C0�u�_ۀ�}$&�  C0zu��4c�}G��  C0�T N���}�"�  C0��m0�} ȡ   C0��4̼^�|�	   C0�>��a��|�I�@  C0�O�L��|݊`  C0�@��*?�|�ʙ�  C1���Ԋ�|��  C19�m���|�K��  C1^����v�|���  C1��_���|�̒   C1�D�i��|�   C1­����|�M�@  C1�),����|�`  C1�����|sΊ�  C2��g���|h�  C2,�}�|\O��  C2?+G6��|P��  C2Rchr$_�|DЃ   C2c��v.T�|9   C2r�skj��|-Q@  C2�����`�|!��`  C2���C��|�{�  C2�~x�t�|
��  C2�~L�{�Sw�  C2��7O,�{���  C2�|�/A�{��t   C2�\^�@�{��   C3
j��{�Up@  C3+�|X.L�{Õ�`  C3K�<�o�{��l�  C3i���i�{��  C3�$W'���{�Wh�  C3�kȫ��{����  C3�;���.�{��e   C3�vn6��{}�   C3���O�{qYa@  C3�!�����{e��`  C3�<��b�{Y�]�  C4OڧT��{N۠  C4�QDJ:�{B[Y�  C4�J"m[�{6���  C4�x����{*�V   C4$���{�   C4,�Q�2x�{]R@  C46!Av���{��`  C4H��k��z��N�  C4c�"Б�z�̠  C4���W&�z�_J�  C4��CG��z؟��  C4�i({BN�z��G   C4ϒ����z� �   C4����z�aC@  C4�pTw#0�z���`  C3O�˿}e�z��?�  C2�?�cy�z�"��  C1��I܃��z�c;�  C1Y�-����zz���  C1�����zn�8   C0Π�����zc$�   C0��F1H�zWe4@  C0��F@���zK��`  C0i���d�z?�0�  C0UQ�z��z4&��  C0TOE��z(g,�  C0n!bB:Z�z���  C0����=(�z�)   C0���+7��z(�   C0��_���y�i%@  C1��1��y���`  C1,��U��y��!�  C1P?���y�*��  C1p]��m �y�k�  C1�w}���y����  C1���¾�y��   C1�
.W���y�,�   C1ٳna�F�y�m@  C1��8�S�y���`  C2��]��y���  C2WY!s��yx.��  C2$y= �P�ylo�  C22��"�-�y`���  C2?w`�  �yT�   C2K6<c��yI0�   C2W)�֮��y=q@  C2f�2=��y1��`  C2}���	��y%��  C2���Q&�y2��  C2��b-�{�yr��  C2ḋeoM�y�}�  C3<{����x���   C3 ��'S/�x�4z   C3=6����x�t�@  C3V�~��"�xӵv`  C3n+���x���  C3��-#T�x�6r�  C3�9��$�x�v��  C3� �!��x��n�  C3�s����x���   C3��$�x�8k   C3��LD�c�x�x�@  C3�5�$�{�xu�g`  C3�y�P(�xi��  C3�@��]�x^:c�  C3��aٰ�xRz��  C3� 9�1��xF�_�  C4kD�<�x:��   C4,g���x/<\   C4I������x#|�@  C4g��F��x�X`  C4�������x�ր  C4�c��P��x >T�  C4�V��*�w�~��  C4��}q��w�P�  C4k9����w���   C3��hO�w�@M   C2D,����wŀ�@  C1�n�{F�w��I`  C1!ʀz���w�ǀ  C0�Q@�S��w�BE�  C0�4,�&��w����  C0tE� �w��A�  C0S�iPd �w�   C09�����wsD>   C0%��r^��wg��@  C0�ٶ��w[�:`  C0+H�ō��wP��  C0X�A g=�wDF6�  C0�Ǭ�n��w8���  C0�vHz(��w,�2�  C0���~�H�w!�   C0����I�wH/   C1#�#i��w	��@  C1=�Q)��v��+`  C1Z��_md�v�	��  C1uZ��1�v�J'�  C1�� 
O��vڊ��  C1��)����v��#�  C1�C�#�8�v��   C1ʄ@����v�L    C1�Ԑ�ɿ�v���@  C1�1�~�v��`  C1�3�*���v���  C2ZY���v�N�  C2�A`���v|���  C2�sɯI�vp��  C25],eAU�ve�   C2U�`���vYP   C2w��8��vM��@  C2�/�l�#�vA�`  C2��nKx�v6��  C2֍�� �v*R	�  C2�eJ7(��v���  C3Jy��#�v��  C3#� rS`�v�   C38��Kz�u�T   C3L�fk��u@  C3`eW�"�u���`  C3u��F�u�|�  C3��:4��u�U��  C3�CI8��u��x�  C3�U��\��u����  C3�OU����u�u   C3���m��u�W�   C3�YQ��^�u��q@  C3�E���Q�u���`  C3�~�c��uzm�  C3�a���unY�  C4�"��ub�i�  C4%�)dP�uV���  C4B �Ҽ��uKf   C4\#�뭴�u?[�   C4s�y���u3�b@  C4�2Ҳ��u'��`  C4������u^�  C3��Q���u]ܠ  C2��AVl�u�Z�  C1����.w�t����  C1hĴ���t�W   C0�:l
��t�_�   C0�e�Jx�tՠS@  C0{E�F�x�t���`  C0R�wHxP�t�!O�  C04U����t�a͠  C0S��#��t��K�  C0�ó��t����  C0 xψ�X�t�#H   C0F<���t�c�   C0y�,��tw�D@  C0�$�e�tk��`  C01=�N���t`%@�  C0[L���tTe��  C0����tH�<�  C0����'b�t<��  C0����b�t1'9   C0�$��F��t%g�   C0�%���t�5@  C1	��#�t�`  C1~қq1�t)1�  C11���G�s�i��  C1A�ѻd��s�-�  C1P������s���  C1^�|��E�s�+*   C1k��|}g�s�k�   C1x\�!G��s��&@  C1�%����s��`  C1�iE��g�s�-"�  C1��[2��s�m��  C1�:�lj��s���  C2�9�*�s���  C2*K7DC�su/   C2I����]�sio�   C2h8|}_o�s]�@  C2��}���sQ�`  C2�[� ���sF1�  C2��Z���s:q��  C2Ʌ.U�s.��  C2��{���s"��  C2��(+J�s3   C2�iS����ss�   C3Q(�#��r��@  C3ތ��r��`  C3�Ľθ�r�5�  C3%���Au�r�u��  C3/,���3�rж �  C3;Lֿ���r��~�  C3P)&ּ��r�6�   C3m��t�r�w{   C3�^�H�r���@  C3���o'�r��w`  C3�D�~��r�8��  C3�1J�&�r~ys�  C3�p����rr���  C4]��DC�rf�o�  C4+~�֊�r[:�   C4=1�wX�rO{l   C4M�<]�rC��@  C4a���j�r7�h`  C4u�6�Y�r,<�  C4}��H���r }d�  C4�\�ԡB�r���  C4�vY�% �r�`�  C4����c�q�>�   C4��mHV�q�]   C4�	�f��q��@  C4��(�OI�q� Y`  C4�4s`_~�q�@׀  C3�6�+���qU�  C2�������q����  C1ޔ� ��q�Q�  C1;������q�B�   C0ӚzF�M�q��N   C0��+j@h�q���@  C0W7���&�q|J`  C0*�&+tN�qpDȀ  C0�6w���qd�F�  C/�^$3�]�qX���  C/�4����qMB�  C/���tc��qAF�   C/y�/�c��q5�?   C/K�w�Q��q)ǽ@  C/B)h^�	�q;`  C/u�����qH��  C/�����q�7�  C/��DǪ�p�ɵ�  C/�*z9��p�
3�  C0������p�J�   C0:sb����p׋0   C0a#���|�p�ˮ@  C0�s�����p�,`  C0��Wp��p�L��  C0��HdM��p��(�  C0�G��p�ͦ�  C1A��y�p�$�  C1-���p�N�   C1H���u��py�!   C1a�hAJ��pmϟ@  C1v{#�7�pb`  C1��H�^��pVP��  C1��F0A��pJ��  C1�9�Ҕf�p>ї�  C1�W]�sC�p3�  C1�MX���p'R�   C1��a��t�p�   C1��$L�pӐ@  C1���yLX�p`  C2P�j��o�   C21U�X�o�*@  C2T�w�k�o���  C2vְq�K�o�,�  C2�4�H�x�o��
   C2�u&��o{.@  C2Ԋ<�]�oc��  C2�tY-�,�oL/��  C3fe~�o4��   C3�s����o1�@  C3/�2fl�o��  C3A!���n�3��  C3P��{�nִ�   C3^���=�n�5�@  C3j�����n���  C3uWY!���n�7��  C3~�h��,�nx��   C3�	�f�I�na9�@  C3�hy�`�nI�Հ  C3��rR��n2;��  C3��w"��n��   C3ˠ4��
�n=�@  C3�R��ym�m�ƀ  C4��#��m�?��  C4$��9���m���   C4>��d��m�A�@  C4T��v� �m�·�  C4aF-T)z�mvC��  C3t+a���m^İ   C2���sq�mGE�@  C1��l'�+�m/ƨ�  C13��z���mG��  C0�ĕ�E��m ȡ   C0q�hF��l�I�@  C06������l�ʙ�  C0@=$F��l�K��  C/�����l�̒   C/��`���l�M�@  C/P�A.^��lsΊ�  C/,j�@3(�l\O��  C/&�l���lDЃ   C/6�����l-Q@  C/N����/�l�{�  C/Z��n���k�Sw�  C/X�����k��t   C/QJ�����k�Up@  C/E�UJ��k��l�  C/4���1��k�Wh�  C/*k?q\��k��e   C/K�:�N��kqYa@  C/�2�<M��kY�]�  C/�B�ak)�kB[Y�  C0>ȷ�J�k*�V   C0��t�y�k]R@  C02r�����j��N�  C0E�����j�_J�  C0T�����j��G   C0f'd 9��j�aC@  C0wa(=���j��?�  C0��`��j�c;�  C0�q���8�jn�8   C0�Ў�G�jWe4@  C0�hTN���j?�0�  C1-iy��j(g,�  C1&����G�j�)   C1HY��2��i�i%@  C1l{a�)��i��!�  C1���a{/�i�k�  C1��"KtD�i��   C1��8B��i�m@  C1�_����i���  C1�z�?;:�ilo�  C1�U�P���iT�   C2~r����i=q@  C2ŀY���i%��  C2)�Ǘ�w�ir��  C26->�u��h���   C2A�#�a��h�t�@  C2M?�%�J�h���  C2[�1�ִ�h�v��  C2r�J���h���   C2��H6t6�h�x�@  C2�nh<�hi��  C2ո�~R�hRz��  C2�N���h:��   C3���q�h#|�@  C3.�;|�h�ր  C3F����g�~��  C3^0Yhf#�g���   C3r�.��gŀ�@  C3��ro��g�ǀ  C3�ǒ�Y��g����  C3���Lt��g�   C3������gg��@  C3�r\��gP��  C3ǈ�)&�g8���  C3���b��g!�   C3ץ:��q�g	��@  C3ߙE�d�f�	��  C3�G=���fڊ��  C3�gէs�f��   C4�q����f���@  C45�#k��f���  C4PR9p��f|���  C4B�����fe�   C3�L���fM��@  C25��ԃ�f6��  C1hy�*��f���  C0��~�f�   C0~qvJǖ�e@  C0<rA�W��e�|�  C0	���e��x�  C/��(����e�u   C/|ڇ�/�e��q@  C/XS�@�H�ezm�  C//l�����eb�i�  C/
�/u��eKf   C.��?���e3�b@  C.��A�/�e^�  C.��}^$I�e�Z�  C.��~"��d�W   C/Q�����dՠS@  C/�����~�d�!O�  C/������d��K�  C0�!���d�#H   C0C�m���dw�D@  C0d��fq�d`%@�  C0�_����dH�<�  C0�;_}�V�d1'9   C0��I]��d�5@  C0ӺN��d)1�  C0�Lā��c�-�  C0��	̓ �c�+*   C1������c��&@  C1#揀��c�-"�  C13X�6���c���  C1ATojd�cu/   C1Nz���*�c]�@  C1\CqR���cF1�  C1oM��Y��c.��  C1�����n�c3   C1��G`m0�b��@  C1������b�5�  C1�)`���bж �  C2�|҇	�b�6�   C2/rʡר�b���@  C2M��7�b�8��  C2h��1�br���  C2�S�\��b[:�   C2���ڼA�bC��@  C2�!G��b,<�  C2��u��b���  C2���!̉�a�>�   C2�}�����a��@  C2�D���a�@׀  C2� �2_�a����  C3Ǽc���a�B�   C3��e�0�a���@  C3o�����apDȀ  C3(Bu*���aX���  C3=��m���aAF�   C3[�0���a)ǽ@  C3z�C�$�aH��  C3�B��a�`�ɵ�  C3�;���`�J�   C3�d�SP��`�ˮ@  C3�(s��`�L��  C3����@�`�ͦ�  C3�6C�6.�`�N�   C2�F�����`mϟ@  C1���D��`VP��  C1M*�"�Y�`>ї�  C0�o*�E�`'R�   C0i	!'��`Ӑ@  C0+�᳗�_�   C/��q(0��_���  C/�|����_��
   C/}�ۺ�_c��  C/U��f=�_4��   C/C�ۈ�J�_��  C/IypO�^ִ�   C/|*�I\�^���  C/����j�^x��   C0��UT#�^I�Հ  C0?7݉�Y�^��   C0d*����]�ƀ  C0����Q��]���   C0�0���]�·�  C0�P�%���]^İ   C0�T���]/ƨ�  C0��h1jn�] ȡ   C1��(7�\�ʙ�  C1#iQc�\�̒   C16N;<���\sΊ�  C1GS�����\DЃ   C1V�6`��\�{�  C1dV����[��t   C1q�3!L�[��l�  C1~
7?��[��e   C1�y�L�[Y�]�  C1��ݵ��[*�V   C1ȳV�[��Z��N�  C1�=-�	�Z��G   C2��8;�Z��?�  C20Ae]4�Zn�8   C2P/b���Z?�0�  C2m��~&�Z�)   C2������Y��!�  C2��(�]��Y��   C2���`��Y���  C2���s��YT�   C2�v�����Y%��  C2�����<�X���   C2�}zcS��X���  C3� ��X���   C3I��� �Xi��  C3H�n��X:��   C3&2�˨3�X�ր  C3/w�n^��W���   C3;����7�W�ǀ  C3P?8�G��W�   C3m���+��WP��  C3����#p�W!�   C3�1̚��V�	��  C3ȁ^�p��V��   C3��f��V���  C3��<�M1�Ve�   C4�S��r�V6��  C4#��O�V�   C3���_��U�|�  C2�d)�ha�U�u   C1�֬^���Uzm�  C1LPL6t��UKf   C0�q�4bT�U^�  C0d�HlH�T�W   C05/~��T�!O�  C/��n3�8�T�#H   C/dA\��T`%@�  C/"���^y�T1'9   C.酡����T)1�  C.��箨��S�+*   C.��T��S�-"�  C.�_H��C�Su/   C.�¿��SF1�  C.�����S3   C.�������R�5�  C.؜i��4�R�6�   C.�e�Mo��R�8��  C.�����R[:�   C.��O6�c�R,<�  C.�,���Q�>�   C/8g_I���Q�@׀  C/s��Ha�Q�B�   C/��2w�b�QpDȀ  C/��^�E�QAF�   C/���H���QH��  C0��N��P�J�   C03�&m�P�L��  C0&�/���P�N�   C06f���PVP��  C0R�r����P'R�   C0t�f�TW�O�   C0�O�����O��
   C0�������O4��   C0�u��4�Nִ�   C1�Yy���Nx��   C1&���N��   C1D�G��l�M���   C1a-\+t��M^İ   C1{����M ȡ   C1����f��L�̒   C1�v>�u��LDЃ   C1�+#,:��K��t   C1�>�o�H�K��e   C1�\r�lo