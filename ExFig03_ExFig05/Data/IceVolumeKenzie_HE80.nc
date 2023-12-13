CDF  ?   
      time       lon       lat          	   CDI       @Climate Data Interface version 2.0.3 (https://mpimet.mpg.de/cdi)   Conventions       CF-1.5     source        #PISM development no-version-control    command      / /home/m/m300792/ModelCodes/mPISM_build/install_03-10-22/bin/pismr -i pism_-065000/pism_-065000.nc -stress_balance ssa+sia -y 100 -ys -65000 -o pism_-064900/pism_-064900.nc -surface_given_file forcing/smb_-064900.nc -surface given -ocean_th_file forcing/oce_-064900.nc -ocean th -topg_file pism_-065000/pism_-065000_new_topg.nc -bed_def prescribed -uplift_file pism_-065000/pism_-065000_dbdt.nc -extra_times 10 -extra_vars thk,usurf,bheatflx,velsurf_mag,velbase_mag,tillwat,bmelt,velbase,flux_divergence,velsurf,velbar,wvelsurf,wvelbase,climatic_mass_balance_cumulative,potential_climatic_mass_balance_cumulative,tempbase,tempsurf,mask,temppabase,tempicethk_basal,topg,nonneg_flux_cumulative,grounded_basal_flux_cumulative,floating_basal_flux_cumulative,discharge_flux_cumulative,Href,taud_mag,lat,lon,calving_mask,sea_level -extra_file pism_-064900/pism_-064900_extra -extra_split -ts_times yearly -ts_file pism_-064900/pism_-064900_ts.nc -topg_to_phi 15.0,30.0,-300.0,400.0 -ssafd_ksp_rtol 1e-7 -enthalpy_temperate_conductivity_ratio 0.001 -pseudo_plastic -pseudo_plastic_q .25 -pseudo_plastic_uthreshold 70 -till_effective_fraction_overburden 0.04 -tauc_slippery_grounding_lines -sia_e 8 -ssa_e 1 -ssa_upwind_factor 0.7 -calving eigen_calving,thickness_calving -thickness_calving_threshold 100 -part_redist -part_grid -cfbc -kill_icebergs -eigen_calving_K 1e16 -verbose 2 -sliding_scale_factor_reduces_tauc 1 -th_heat_coefficient 4e-5 -th_salinity_coefficient 4e-8 -o_format netcdf3 -sediment_file till_mask.nc -no_divide_tillphi -periodicity none -ocean_th_period 1 -o_order yxz
    proj4         j+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs     history_of_appended_files        �Fri Jul  8 17:44:35 2022: Appended file outFinal2 had following "history" attribute:
Fri Jul  8 17:44:34 2022: ncap2 -O -s lon_bnds=float(lon_bnds) outFinal1 outFinal2
Fri Jul  8 17:44:33 2022: ncap2 -O -s lat_bnds=float(lat_bnds) outFinal outFinal1
Fri Jul 08 17:44:32 2022: cdo -add thkNew thkOld outFinal
Fri Jul 08 17:44:30 2022: cdo -mul -selvar,thk pism_-064900/pism_-064900.nc MaskLaurentide thkNew
     NCO       _netCDF Operators version 5.0.6 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)    history       |Thu Sep 15 09:42:09 2022: cdo -mergetime File000001.nc File000002.nc File000003.nc File000004.nc File000005.nc File000006.nc File000007.nc File000008.nc File000009.nc File000010.nc File000011.nc File000012.nc File000013.nc File000014.nc File000015.nc File000016.nc File000017.nc File000018.nc File000019.nc File000020.nc File000021.nc File000022.nc File000023.nc File000024.nc File000025.nc File000026.nc File000027.nc File000028.nc File000029.nc File000030.nc File000031.nc File000032.nc File000033.nc File000034.nc File000035.nc File000036.nc File000037.nc File000038.nc File000039.nc File000040.nc File000041.nc File000042.nc File000043.nc File000044.nc File000045.nc File000046.nc File000047.nc File000048.nc File000049.nc File000050.nc File000051.nc File000052.nc File000053.nc File000054.nc File000055.nc File000056.nc File000057.nc File000058.nc File000059.nc File000060.nc File000061.nc File000062.nc File000063.nc File000064.nc File000065.nc File000066.nc File000067.nc File000068.nc File000069.nc File000070.nc File000071.nc File000072.nc File000073.nc File000074.nc File000075.nc File000076.nc File000077.nc File000078.nc File000079.nc File000080.nc File000081.nc File000082.nc File000083.nc File000084.nc File000085.nc File000086.nc File000087.nc File000088.nc File000089.nc File000090.nc File000091.nc File000092.nc File000093.nc File000094.nc File000095.nc File000096.nc File000097.nc File000098.nc File000099.nc File000100.nc File000101.nc File000102.nc File000103.nc File000104.nc File000105.nc File000106.nc File000107.nc File000108.nc File000109.nc File000110.nc File000111.nc File000112.nc File000113.nc File000114.nc File000115.nc File000116.nc File000117.nc File000118.nc File000119.nc File000120.nc File000121.nc File000122.nc File000123.nc File000124.nc File000125.nc File000126.nc File000127.nc File000128.nc File000129.nc File000130.nc File000131.nc File000132.nc File000133.nc File000134.nc File000135.nc File000136.nc File000137.nc File000138.nc File000139.nc File000140.nc File000141.nc File000142.nc File000143.nc File000144.nc File000145.nc File000146.nc File000147.nc File000148.nc File000149.nc File000150.nc File000151.nc File000152.nc File000153.nc File000154.nc File000155.nc File000156.nc File000157.nc File000158.nc File000159.nc File000160.nc File000161.nc File000162.nc File000163.nc File000164.nc File000165.nc File000166.nc File000167.nc File000168.nc File000169.nc File000170.nc File000171.nc File000172.nc File000173.nc File000174.nc File000175.nc File000176.nc File000177.nc File000178.nc File000179.nc File000180.nc File000181.nc File000182.nc File000183.nc File000184.nc File000185.nc File000186.nc File000187.nc File000188.nc File000189.nc File000190.nc File000191.nc File000192.nc File000193.nc File000194.nc File000195.nc File000196.nc File000197.nc File000198.nc File000199.nc File000200.nc File000201.nc File000202.nc File000203.nc File000204.nc File000205.nc File000206.nc File000207.nc File000208.nc File000209.nc File000210.nc File000211.nc File000212.nc File000213.nc File000214.nc File000215.nc File000216.nc File000217.nc File000218.nc File000219.nc File000220.nc File000221.nc File000222.nc File000223.nc File000224.nc File000225.nc File000226.nc File000227.nc File000228.nc File000229.nc File000230.nc File000231.nc File000232.nc File000233.nc File000234.nc File000235.nc File000236.nc File000237.nc File000238.nc File000239.nc File000240.nc File000241.nc File000242.nc File000243.nc File000244.nc File000245.nc File000246.nc File000247.nc File000248.nc File000249.nc File000250.nc File000251.nc File000252.nc File000253.nc File000254.nc File000255.nc File000256.nc File000257.nc File000258.nc File000259.nc File000260.nc File000261.nc File000262.nc File000263.nc File000264.nc File000265.nc File000266.nc File000267.nc File000268.nc File000269.nc File000270.nc File000271.nc File000272.nc File000273.nc File000274.nc File000275.nc File000276.nc File000277.nc File000278.nc File000279.nc File000280.nc File000281.nc File000282.nc File000283.nc File000284.nc File000285.nc File000286.nc File000287.nc File000288.nc File000289.nc File000290.nc File000291.nc File000292.nc File000293.nc File000294.nc File000295.nc File000296.nc File000297.nc File000298.nc File000299.nc File000300.nc File000301.nc File000302.nc File000303.nc File000304.nc File000305.nc File000306.nc File000307.nc File000308.nc File000309.nc File000310.nc File000311.nc File000312.nc File000313.nc File000314.nc File000315.nc File000316.nc File000317.nc File000318.nc File000319.nc File000320.nc File000321.nc File000322.nc File000323.nc File000324.nc File000325.nc File000326.nc File000327.nc File000328.nc File000329.nc File000330.nc File000331.nc File000332.nc File000333.nc File000334.nc File000335.nc File000336.nc File000337.nc File000338.nc File000339.nc File000340.nc File000341.nc File000342.nc File000343.nc File000344.nc File000345.nc File000346.nc File000347.nc File000348.nc File000349.nc File000350.nc File000351.nc File000352.nc File000353.nc File000354.nc File000355.nc File000356.nc File000357.nc File000358.nc File000359.nc File000360.nc File000361.nc File000362.nc File000363.nc File000364.nc File000365.nc File000366.nc File000367.nc File000368.nc File000369.nc File000370.nc File000371.nc File000372.nc File000373.nc File000374.nc File000375.nc File000376.nc File000377.nc File000378.nc File000379.nc File000380.nc File000381.nc File000382.nc File000383.nc File000384.nc File000385.nc File000386.nc File000387.nc File000388.nc File000389.nc File000390.nc File000391.nc File000392.nc File000393.nc File000394.nc File000395.nc File000396.nc File000397.nc File000398.nc File000399.nc File000400.nc File000401.nc File000402.nc File000403.nc File000404.nc File000405.nc File000406.nc File000407.nc File000408.nc File000409.nc File000410.nc File000411.nc File000412.nc File000413.nc File000414.nc File000415.nc File000416.nc File000417.nc File000418.nc File000419.nc File000420.nc File000421.nc File000422.nc File000423.nc File000424.nc File000425.nc File000426.nc File000427.nc File000428.nc File000429.nc File000430.nc File000431.nc File000432.nc File000433.nc File000434.nc File000435.nc File000436.nc File000437.nc File000438.nc File000439.nc File000440.nc File000441.nc File000442.nc File000443.nc File000444.nc File000445.nc File000446.nc File000447.nc File000448.nc File000449.nc File000450.nc File000451.nc File000452.nc File000453.nc File000454.nc File000455.nc File000456.nc File000457.nc File000458.nc File000459.nc File000460.nc File000461.nc File000462.nc File000463.nc File000464.nc File000465.nc File000466.nc File000467.nc File000468.nc File000469.nc File000470.nc File000471.nc File000472.nc File000473.nc File000474.nc File000475.nc File000476.nc File000477.nc File000478.nc File000479.nc File000480.nc File000481.nc File000482.nc File000483.nc File000484.nc File000485.nc File000486.nc File000487.nc File000488.nc File000489.nc File000490.nc File000491.nc File000492.nc File000493.nc File000494.nc File000495.nc File000496.nc File000497.nc File000498.nc File000499.nc File000500.nc File000501.nc File000502.nc File000503.nc File000504.nc File000505.nc File000506.nc File000507.nc File000508.nc File000509.nc File000510.nc File000511.nc File000512.nc File000513.nc File000514.nc File000515.nc File000516.nc File000517.nc File000518.nc File000519.nc File000520.nc File000521.nc File000522.nc File000523.nc File000524.nc File000525.nc File000526.nc File000527.nc File000528.nc File000529.nc File000530.nc File000531.nc File000532.nc File000533.nc File000534.nc File000535.nc File000536.nc File000537.nc File000538.nc File000539.nc File000540.nc File000541.nc File000542.nc File000543.nc File000544.nc File000545.nc File000546.nc File000547.nc File000548.nc File000549.nc File000550.nc File000551.nc File000552.nc File000553.nc File000554.nc File000555.nc File000556.nc File000557.nc File000558.nc File000559.nc File000560.nc File000561.nc File000562.nc File000563.nc File000564.nc File000565.nc File000566.nc File000567.nc File000568.nc File000569.nc File000570.nc File000571.nc File000572.nc File000573.nc File000574.nc File000575.nc IceVolumeKenzie.nc
Thu Sep 15 09:38:14 2022: cdo -fldsum -mulc,10000 -mulc,10000 -selindexbox,211,427,600,415 -selvar,thk /work/ba0989/m300792/HE_Runs/ExperimentsComposite//HE80//pism_-064900/pism_-064900.nc Tmp/File000001.nc   CDO       @Climate Data Operators version 2.0.3 (https://mpimet.mpg.de/cdo)         time                standard_name         time   	long_name         time   units         seconds since 1-1-1    calendar      365_day    axis      T               -�   lon                standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X               -�   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y               -�   thk                       standard_name         land_ice_thickness     	long_name         land ice thickness     units         m      pism_intent       model_state             -�                �}Ȁ@�  C-T�I�T�}���   C-�(z[�}�=   C-��R�f�}�A�@  C.��&"}�}��9`  C.!�v�k��}�·�  C.:&����}�5�  C.X�B+r��}vC��  C.m0M��#�}j�1�  C.x�A:�H�}^İ   C.����3��}S.   C.��}�^�}GE�@  C/;&뎻E�};�*`  C/�nV�u�}/ƨ�  C0��X�}$&�  C0.�����}G��  C0 �����}�"�  C0�QL*��} ȡ   C/�%|c"�|�	   C/�@,l<��|�I�@  C/��@���|݊`  C/�F���|�ʙ�  C02�p�+��|��  C0lp�L�|�K��  C0��w�O�|���  C0�&}����|�̒   C0�m:���|�   C0٣:㈃�|�M�@  C0�����|�`  C0�ί?%�|sΊ�  C0�w *�\�|h�  C1L!�S��|\O��  C1@�+	/��|P��  C1s�y�07�|DЃ   C1�O?d���|9   C1���7��|-Q@  C1����d�|!��`  C1��B���|�{�  C1n�`�<�|
��  C1Q�[:�{�Sw�  C0�q����{���  C/��[��{��t   C.�b���{��   C.B{b�G��{�Up@  C-�&���{Õ�`  C-=(�+��{��l�  C,ꮾ���{��  C,�� �v4�{�Wh�  C,���BȤ�{����  C,qM����{��e   C,f�|oh��{}�   C,�RxB��{qYa@  C-"I1=.�{e��`  C-�::�U��{Y�]�  C-����W��{N۠  C.(�t��V�{B[Y�  C.T�����{6���  C.wuB�ּ�{*�V   C.5�(y���{�   C-�v��=W�{]R@  C-�.y]��{��`  C-���:���z��N�  C.C�����z�̠  C.����z�z�_J�  C.���{���z؟��  C/Zp��;�z��G   C/u�� g�z� �   C/@�q?�z�aC@  C/bt�K�z���`  C/�l�z��?�  C/�F0;�z�"��  C/S�����z�c;�  C/�P*^��zz���  C0"��i�zn�8   C0U�	3�r�zc$�   C0�VX]��zWe4@  C0����8�zK��`  C0���C4��z?�0�  C0��� ��z4&��  C0����~�z(g,�  C0���jކ�z���  C0�#k'��z�)   C0���é�z(�   C0�|���y�i%@  C1��=�!�y���`  C11F�����y��!�  C1>��N�y�*��  C1@�.WZr�y�k�  C0��� sT�y����  C/���TO-�y��   C.�I� b��y�,�   C.c=�<���y�m@  C.1��Y��y���`  C.;����y���  C-Օ�A��yx.��  C-�A�Kd��ylo�  C-|��GQ��y`���  C- K��]%�yT�   C,� �VI��yI0�   C,]��H�y=q@  C,>q��p��y1��`  C,l�[��^�y%��  C,ک뭴��y2��  C-D���X�yr��  C-��@���y�}�  C-�Ֆ&��x���   C. �p����x�4z   C.G�eA<�x�t�@  C.`Q�*P=�xӵv`  C.fI��AN�x���  C.i*٨��x�6r�  C.1T�\z�x�v��  C.d� "���x��n�  C.�x��w�x���   C.ש}�F��x�8k   C/P)#� �x�x�@  C/p�n��xu�g`  C/'������xi��  C/*�1Ip�x^:c�  C/0*ġ���xRz��  C/:.���xF�_�  C/}Q�Q���x:��   C/� �V�x/<\   C061�?���x#|�@  C0y\�ߨ�x�X`  C0�X��>��x�ր  C0�i��$�x >T�  C0�XIھ�w�~��  C0�S۾���w�P�  C0�^P]���w���   C0��]�z;�w�@M   C0��6UN��wŀ�@  C0αQ�M�w��I`  C0�k� k�w�ǀ  C1�s�'�w�BE�  C1&�P���w����  C1+�<'��w��A�  C1!}�
��w�   C15�@5��wsD>   C1�\k��wg��@  C0��ci���w[�:`  C/&/O��2�wP��  C.���_@��wDF6�  C.6,-�N�w8���  C-�^FՊ��w,�2�  C-~T�2QJ�w!�   C-�o��v�wH/   C,����u��w	��@  C,x¸����v��+`  C,f����v�	��  C,W[{,=�v�J'�  C,�0hC���vڊ��  C,�ᒪ���v��#�  C-e����v��   C-�}��v�L    C.������v���@  C./�i�8��v��`  C.N+	'�v���  C.Y�P�"��v�N�  C.Sgޔ�!�v|���  C.TǴ8�vp��  C.2H���ve�   C.4��d���vYP   C.jN�o��vM��@  C.��O�$)�vA�`  C.ćˌK��v6��  C.ۗ�W�*�v*R	�  C.���D��v���  C/�X��v��  C/�����v�   C/�忛�u�T   C/_�MS��u@  C/��^P�u���`  C09Xr�;7�u�|�  C0| �Ab��u�U��  C0�w�� �u��x�  C0��I���u����  C0�,����u�u   C0�/u�N�u�W�   C0�P&B	��u��q@  C0�D����u���`  C0��{�&n�uzm�  C0�:N���unY�  C0���z��ub�i�  C1��3�N�uV���  C1&��f'�uKf   C1+�}a��u?[�   C1-�w��u3�b@  C1(��:��u'��`  C1:���)�u^�  C1"�9���u]ܠ  C14n��M�u�Z�  C1w�h���t����  C0���
�t�W   C0h�d��t�_�   C/ƙ ���tՠS@  C.!4�X��t���`  C-LMW�m�t�!O�  C,�+荰Q�t�a͠  C,/~�kl�t��K�  C+��`���t����  C,<bX�9�t�#H   C,�ċ-��t�c�   C-���e�tw�D@  C-l�W����tk��`  C-�aR&���t`%@�  C-��w
���tTe��  C-�������tH�<�  C.��ӡ\�t<��  C.�d ��t1'9   C."��$��t%g�   C-�r\1fH�t�5@  C.!|E���t�`  C.`�� �k�t)1�  C.��6��x�s�i��  C.�+�|6��s�-�  C.�4���g�s���  C.�u��{�s�+*   C.������s�k�   C.�tb-,H�s��&@  C.���B3p�s��`  C/4�����s�-"�  C/��x��s�m��  C0�,P'i�s���  C0A`?�l�s���  C0w��R.�su/   C0�՚���sio�   C0��J:�1�s]�@  C0������sQ�`  C0t�`�l�sF1�  C0g]�����s:q��  C0z�/�Ą�s.��  C0�������s"��  C0�}���s3   C0�&�{�ss�   C1\0{��r��@  C18Ln��r��`  C1n:�8�r�5�  C1R����r�u��  C1��e1~�rж �  C1�t`%�r��~�  C1(��r=��r�6�   C0i���Q��r�w{   C/�d.u(�r���@  C/mbE|!�r��w`  C.��!_�	�r�8��  C-�v;^E��r~ys�  C-F��eY�rr���  C,�J���rf�o�  C,k�L�$�r[:�   C,?�k��	�rO{l   C,j[����rC��@  C,�4��K��r7�h`  C-(��
��r,<�  C-oZ��P.�r }d�  C-�v\
���r���  C-�襯}��r�`�  C. �x�c�q�>�   C.+�e�7�q�]   C-� !9EX�q��@  C-��)�\��q� Y`  C-��2�(��q�@׀  C-�t�ɀ��qU�  C.��ȉ��q����  C.5������q�Q�  C.^���q�B�   C.~F���&�q��N   C.���}b��q���@  C.�f����q|J`  C.�<�	�qpDȀ  C.�D�'��qd�F�  C.�äD�m�qX���  C/c[���qMB�  C/�ٵo�[�qAF�   C0 ��p�q5�?   C0L������q)ǽ@  C0l�B����q;`  C0Z?�q<�qH��  C0J5�A0}�q�7�  C0@�
2P�p�ɵ�  C09�����p�
3�  C0R����p�J�   C0���-ϯ�p׋0   C0������p�ˮ@  C0���B��p�,`  C0�K9���p�L��  C1�o{��p��(�  C1�����p�ͦ�  C1��jC�p�$�  C1��|c��p�N�   C1��{�py�!   C1$��f�x�pmϟ@  C1=�i+|�pb`  C0���"�pVP��  C0!m2���pJ��  C/������p>ї�  C.hS���p3�  C-vvQ���p'R�   C,�E.tc1�p�   C,9n_V��pӐ@  C+�x��I�p`  C+s�]W��o�   C+�������o�*@  C+ː'b�o���  C,D���D�o�,�  C,e��q���o��
   C,�ϕ�Ԟ�o{.@  C,��͘1��oc��  C,���%?�oL/��  C,�|q���o4��   C,�AW�2��o1�@  C-	���o��  C-~�oS.B�n�3��  C-�Ե]}�nִ�   C.E.�G�n�5�@  C.`7N8d��n���  C.:�C���n�7��  C.,{�B)��nx��   C.��G���na9�@  C.�����nI�Հ  C.���W�n2;��  C.R�;Ʈ��n��   C.�'az�m�n=�@  C/=E%Oq>�m�ƀ  C/�ܛٗ`�m�?��  C0�YUj��m���   C00��j%��m�A�@  C0F�X�m�·�  C0����d�mvC��  C0�ގ��m^İ   C0 ]����mGE�@  C0С۰|�m/ƨ�  C0E}ԥm��mG��  C0i� �(�m ȡ   C0���r���l�I�@  C0���7|�l�ʙ�  C0�͖��l�K��  C0�oxE6j�l�̒   C0������l�M�@  C0�_8DK��lsΊ�  C0�������l\O��  C0I�<w���lDЃ   C/���n&��l-Q@  C/�˰���l�{�  C.���	�k�Sw�  C.X�̖���k��t   C-���7�k�Up@  C-,-�(��k��l�  C,�&����k�Wh�  C,8-u����k��e   C,��"˫�kqYa@  C,&��kY�]�  C,�����`�kB[Y�  C,�:�|���k*�V   C-&>Vd��k]R@  C-V��;��j��N�  C-p-G:d�j�_J�  C-y/�	�j��G   C-{�Ê�-�j�aC@  C-`�loyk�j��?�  C-J��.�j�c;�  C-�(R~��jn�8   C-Y�윹P�jWe4@  C-��� q��j?�0�  C-�A�J�j(g,�  C-� ���i�j�)   C.����T�i�i%@  C.	��Q��i��!�  C.�N�i�k�  C.�X1}y�i��   C.��׶�i�m@  C.Llr=S��i���  C.��,��ilo�  C/4Ϯ�iT�   C/�{C�8/�i=q@  C/��'.���i%��  C0ӿ(Z�ir��  C0��A(��h���   C/�Uңv��h�t�@  C/�e���h���  C/� PY.��h�v��  C/�p'8��h���   C0")n���h�x�@  C0V��'��hi��  C0��Zzb�hRz��  C0��k&���h:��   C0��a����h#|�@  C0��v��h�ր  C0�Q�����g�~��  C0�wu���g���   C0���f�gŀ�@  C0t�A��g�ǀ  C/�L�ѳ��g����  C/����g�   C.�ٳ���gg��@  C.T蕀��gP��  C-�6��k��g8���  C,�Rd�o?�g!�   C,n�,��g	��@  C,���?M�f�	��  C+��x��v�fڊ��  C,�Y��f��   C,r����&�f���@  C,��0��f���  C-@0���f|���  C-�V���<�fe�   C-��pV2��fM��@  C-���:�f6��  C-�ڤ0n�f���  C-��;���f�   C-D�u��e@  C-PÁ�+��e�|�  C-�06#���e��x�  C-���˻�e�u   C.*�?|&r�e��q@  C.^���{��ezm�  C.ykSZ�eb�i�  C.���}�a�eKf   C.����y�e3�b@  C.��#����e^�  C.���o��e�Z�  C.�
h`'��d�W   C/[��<�dՠS@  C/�XO���d�!O�  C0&�!�~�d��K�  C0\�}��3�d�#H   C0[	���9�dw�D@  C0Y��t��d`%@�  C0X6����dH�<�  C0U����d1'9   C0Rɗ� B�d�5@  C0kf�6P��d)1�  C0�j����c�-�  C0��#�c�+*   C0�]=�c��&@  C0�^?�c�-"�  C0`a�d�c���  C/7eы!��cu/   C.EX�B���c]�@  C-�/\���cF1�  C-L�d��c.��  C,�p?f��c3   C,�z����b��@  C,���J�b�5�  C,�y����bж �  C,�
��
�b�6�   C,t��W�b���@  C,"ŻB�A�b�8��  C+��TU"�br���  C+�2H0E��b[:�   C+��w����bC��@  C,&+�����b,<�  C,�A<A��b���  C-
-��a�>�   C-h)��E�a��@  C-���ol�a�@׀  C-��o��a����  C.�$~k��a�B�   C.*�s�)4�a���@  C-��ye���apDȀ  C-�PL(o��aX���  C-��O���aAF�   C.�zv�c�a)ǽ@  C.Y$�;V�aH��  C.�rP� ��`�ɵ�  C.�B����`�J�   C.�՝����`�ˮ@  C.�e&n�`�L��  C.�k���z�`�ͦ�  C.�}jt3�`�N�   C/y:�iH�`mϟ@  C/l�~�=��`VP��  C0�gE��`>ї�  C0[�0r4��`'R�   C0��Yǌ�`Ӑ@  C0�f@F���_�   C0���!��_���  C0���T��_��
   C0�!4b���_c��  C0.zK7���_4��   C/>щ��_��  C.N�f!�^ִ�   C-���g�^���  C-MHl�6�^x��   C,�^7�)x�^I�Հ  C,��]���^��   C,L���-`�]�ƀ  C+�{�|}K�]���   C+��f��E�]�·�  C+[ѤR5�]^İ   C+:`"e��]/ƨ�  C+eg�25��] ȡ   C+ƯN>���\�ʙ�  C,#q=�M��\�̒   C,l?6���\sΊ�  C,�aq����\DЃ   C,�-!.���\�{�  C,ί\��=�[��t   C,�ݗ&��[��l�  C,M��ځ��[��e   C,&���{�[Y�]�  C,ZAu�S��[*�V   C,�z�Z�Z��N�  C-Yp�!#�Z��G   C-����T�Z��?�  C.!=;��*�Zn�8   C.]��2�Z?�0�  C.�l&L{��Z�)   C.�]�N��Y��!�  C.�|�����Y��   C.�BZ*6��Y���  C.�  V;��YT�   C.�7��J��Y%��  C/ m����X���   C/5�M�P[�X���  C/f�ɀ��X���   C/��BR�H�Xi��  C/�7��#&�X:��   C/�[��-�X�ր  C/����v�W���   C/���V�W�ǀ  C0��[X��W�   C0aDSE���WP��  C0��O�hA�W!�   C0ֲ Z�V�	��  C0���I[<�V��   C/�t՛�x�V���  C.UL���Ve�   C-���b�l�V6��  C,�vG�D�V�   C,`�a�j�U�|�  C,!����z�U�u   C,�/8�Uzm�  C+�@�����UKf   C+��j���U^�  C+��ԧD�T�W   C+�
T����T�!O�  C+��c���T�#H   C+f=�K�f�T`%@�  C+J2��!Z�T1'9   C+[�����T)1�  C+��J,�S�+*   C,^ ����S�-"�  C,�[D����Su/   C,�U��G�SF1�  C-�v����S3   C-S��R��R�5�  C-o��	���R�6�   C-A��b���R�8��  C-D©��R[:�   C,쥡d���R,<�  C-d����Q�>�   C-kn��.�Q�@׀  C-�g�Y��Q�B�   C-��MK:��QpDȀ  C..�,b��QAF�   C.C�z��QH��  C.L�ٲ���P�J�   C.QCf����P�L��  C.^�+���P�N�   C.w�#��C�PVP��  C.�G�)���P'R�   C/���Vz��O�   C0�#���O��
   C0Pml,)�O4��   C0��9 �Nִ�   C0�Yu��Nx��   C0�x.k���N��   C0w{�����M���   C0m��\��M^İ   C0f�4����M ȡ   C0}��o_�L�̒   C0�)�_���LDЃ   C0��k���K��t   C0�����K��e   C0#~�%P