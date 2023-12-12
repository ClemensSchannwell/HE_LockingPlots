CDF  ?   
      time       lon       lat          	   CDI       @Climate Data Interface version 2.0.3 (https://mpimet.mpg.de/cdi)   Conventions       CF-1.5     source        #PISM development no-version-control    command      / /home/m/m300792/ModelCodes/mPISM_build/install_03-10-22/bin/pismr -i pism_-065000/pism_-065000.nc -stress_balance ssa+sia -y 100 -ys -65000 -o pism_-064900/pism_-064900.nc -surface_given_file forcing/smb_-064900.nc -surface given -ocean_th_file forcing/oce_-064900.nc -ocean th -topg_file pism_-065000/pism_-065000_new_topg.nc -bed_def prescribed -uplift_file pism_-065000/pism_-065000_dbdt.nc -extra_times 10 -extra_vars thk,usurf,bheatflx,velsurf_mag,velbase_mag,tillwat,bmelt,velbase,flux_divergence,velsurf,velbar,wvelsurf,wvelbase,climatic_mass_balance_cumulative,potential_climatic_mass_balance_cumulative,tempbase,tempsurf,mask,temppabase,tempicethk_basal,topg,nonneg_flux_cumulative,grounded_basal_flux_cumulative,floating_basal_flux_cumulative,discharge_flux_cumulative,Href,taud_mag,lat,lon,calving_mask,sea_level -extra_file pism_-064900/pism_-064900_extra -extra_split -ts_times yearly -ts_file pism_-064900/pism_-064900_ts.nc -topg_to_phi 15.0,30.0,-300.0,400.0 -ssafd_ksp_rtol 1e-7 -enthalpy_temperate_conductivity_ratio 0.001 -pseudo_plastic -pseudo_plastic_q .25 -pseudo_plastic_uthreshold 70 -till_effective_fraction_overburden 0.04 -tauc_slippery_grounding_lines -sia_e 8 -ssa_e 1 -ssa_upwind_factor 0.7 -calving eigen_calving,thickness_calving -thickness_calving_threshold 100 -part_redist -part_grid -cfbc -kill_icebergs -eigen_calving_K 1e16 -verbose 2 -sliding_scale_factor_reduces_tauc 1 -th_heat_coefficient 4e-5 -th_salinity_coefficient 4e-8 -o_format netcdf3 -sediment_file till_mask.nc -no_divide_tillphi -periodicity none -ocean_th_period 1 -o_order yxz
    proj4         j+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs     history_of_appended_files        �Tue Apr 26 17:25:11 2022: Appended file outFinal2 had following "history" attribute:
Tue Apr 26 17:25:10 2022: ncap2 -O -s lon_bnds=float(lon_bnds) outFinal1 outFinal2
Tue Apr 26 17:25:10 2022: ncap2 -O -s lat_bnds=float(lat_bnds) outFinal outFinal1
Tue Apr 26 17:25:09 2022: cdo -add thkNew thkOld outFinal
Tue Apr 26 17:25:09 2022: cdo -mul -selvar,thk pism_-064900/pism_-064900.nc MaskLaurentide thkNew
     NCO       _netCDF Operators version 5.0.6 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)    history       lTue Sep 13 22:45:52 2022: cdo -mergetime File000001.nc File000002.nc File000003.nc File000004.nc File000005.nc File000006.nc File000007.nc File000008.nc File000009.nc File000010.nc File000011.nc File000012.nc File000013.nc File000014.nc File000015.nc File000016.nc File000017.nc File000018.nc File000019.nc File000020.nc File000021.nc File000022.nc File000023.nc File000024.nc File000025.nc File000026.nc File000027.nc File000028.nc File000029.nc File000030.nc File000031.nc File000032.nc File000033.nc File000034.nc File000035.nc File000036.nc File000037.nc File000038.nc File000039.nc File000040.nc File000041.nc File000042.nc File000043.nc File000044.nc File000045.nc File000046.nc File000047.nc File000048.nc File000049.nc File000050.nc File000051.nc File000052.nc File000053.nc File000054.nc File000055.nc File000056.nc File000057.nc File000058.nc File000059.nc File000060.nc File000061.nc File000062.nc File000063.nc File000064.nc File000065.nc File000066.nc File000067.nc File000068.nc File000069.nc File000070.nc File000071.nc File000072.nc File000073.nc File000074.nc File000075.nc File000076.nc File000077.nc File000078.nc File000079.nc File000080.nc File000081.nc File000082.nc File000083.nc File000084.nc File000085.nc File000086.nc File000087.nc File000088.nc File000089.nc File000090.nc File000091.nc File000092.nc File000093.nc File000094.nc File000095.nc File000096.nc File000097.nc File000098.nc File000099.nc File000100.nc File000101.nc File000102.nc File000103.nc File000104.nc File000105.nc File000106.nc File000107.nc File000108.nc File000109.nc File000110.nc File000111.nc File000112.nc File000113.nc File000114.nc File000115.nc File000116.nc File000117.nc File000118.nc File000119.nc File000120.nc File000121.nc File000122.nc File000123.nc File000124.nc File000125.nc File000126.nc File000127.nc File000128.nc File000129.nc File000130.nc File000131.nc File000132.nc File000133.nc File000134.nc File000135.nc File000136.nc File000137.nc File000138.nc File000139.nc File000140.nc File000141.nc File000142.nc File000143.nc File000144.nc File000145.nc File000146.nc File000147.nc File000148.nc File000149.nc File000150.nc File000151.nc File000152.nc File000153.nc File000154.nc File000155.nc File000156.nc File000157.nc File000158.nc File000159.nc File000160.nc File000161.nc File000162.nc File000163.nc File000164.nc File000165.nc File000166.nc File000167.nc File000168.nc File000169.nc File000170.nc File000171.nc File000172.nc File000173.nc File000174.nc File000175.nc File000176.nc File000177.nc File000178.nc File000179.nc File000180.nc File000181.nc File000182.nc File000183.nc File000184.nc File000185.nc File000186.nc File000187.nc File000188.nc File000189.nc File000190.nc File000191.nc File000192.nc File000193.nc File000194.nc File000195.nc File000196.nc File000197.nc File000198.nc File000199.nc File000200.nc File000201.nc File000202.nc File000203.nc File000204.nc File000205.nc File000206.nc File000207.nc File000208.nc File000209.nc File000210.nc File000211.nc File000212.nc File000213.nc File000214.nc File000215.nc File000216.nc File000217.nc File000218.nc File000219.nc File000220.nc File000221.nc File000222.nc File000223.nc File000224.nc File000225.nc File000226.nc File000227.nc File000228.nc File000229.nc File000230.nc File000231.nc File000232.nc File000233.nc File000234.nc File000235.nc File000236.nc File000237.nc File000238.nc File000239.nc File000240.nc File000241.nc File000242.nc File000243.nc File000244.nc File000245.nc File000246.nc File000247.nc File000248.nc File000249.nc File000250.nc File000251.nc File000252.nc File000253.nc File000254.nc File000255.nc File000256.nc File000257.nc File000258.nc File000259.nc File000260.nc File000261.nc File000262.nc File000263.nc File000264.nc File000265.nc File000266.nc File000267.nc File000268.nc File000269.nc File000270.nc File000271.nc File000272.nc File000273.nc File000274.nc File000275.nc File000276.nc File000277.nc File000278.nc File000279.nc File000280.nc File000281.nc File000282.nc File000283.nc File000284.nc File000285.nc File000286.nc File000287.nc File000288.nc File000289.nc File000290.nc File000291.nc File000292.nc File000293.nc File000294.nc File000295.nc File000296.nc File000297.nc File000298.nc File000299.nc File000300.nc File000301.nc File000302.nc File000303.nc File000304.nc File000305.nc File000306.nc File000307.nc File000308.nc File000309.nc File000310.nc File000311.nc File000312.nc File000313.nc File000314.nc File000315.nc File000316.nc File000317.nc File000318.nc File000319.nc File000320.nc File000321.nc File000322.nc File000323.nc File000324.nc File000325.nc File000326.nc File000327.nc File000328.nc File000329.nc File000330.nc File000331.nc File000332.nc File000333.nc File000334.nc File000335.nc File000336.nc File000337.nc File000338.nc File000339.nc File000340.nc File000341.nc File000342.nc File000343.nc File000344.nc File000345.nc File000346.nc File000347.nc File000348.nc File000349.nc File000350.nc File000351.nc File000352.nc File000353.nc File000354.nc File000355.nc File000356.nc File000357.nc File000358.nc File000359.nc File000360.nc File000361.nc File000362.nc File000363.nc File000364.nc File000365.nc File000366.nc File000367.nc File000368.nc File000369.nc File000370.nc File000371.nc File000372.nc File000373.nc File000374.nc File000375.nc File000376.nc File000377.nc File000378.nc File000379.nc File000380.nc File000381.nc File000382.nc File000383.nc File000384.nc File000385.nc File000386.nc File000387.nc File000388.nc File000389.nc File000390.nc File000391.nc File000392.nc File000393.nc File000394.nc File000395.nc File000396.nc File000397.nc File000398.nc File000399.nc File000400.nc File000401.nc File000402.nc File000403.nc File000404.nc File000405.nc File000406.nc File000407.nc File000408.nc File000409.nc File000410.nc File000411.nc File000412.nc File000413.nc File000414.nc File000415.nc File000416.nc File000417.nc File000418.nc File000419.nc File000420.nc File000421.nc File000422.nc File000423.nc File000424.nc File000425.nc File000426.nc File000427.nc File000428.nc File000429.nc File000430.nc File000431.nc File000432.nc File000433.nc File000434.nc File000435.nc File000436.nc File000437.nc File000438.nc File000439.nc File000440.nc File000441.nc File000442.nc File000443.nc File000444.nc File000445.nc File000446.nc File000447.nc File000448.nc File000449.nc File000450.nc File000451.nc File000452.nc File000453.nc File000454.nc File000455.nc File000456.nc File000457.nc File000458.nc File000459.nc File000460.nc File000461.nc File000462.nc File000463.nc File000464.nc File000465.nc File000466.nc File000467.nc File000468.nc File000469.nc File000470.nc File000471.nc File000472.nc File000473.nc File000474.nc File000475.nc File000476.nc File000477.nc File000478.nc File000479.nc File000480.nc File000481.nc File000482.nc File000483.nc File000484.nc File000485.nc File000486.nc File000487.nc File000488.nc File000489.nc File000490.nc File000491.nc File000492.nc File000493.nc File000494.nc File000495.nc File000496.nc File000497.nc File000498.nc File000499.nc File000500.nc File000501.nc File000502.nc File000503.nc File000504.nc File000505.nc File000506.nc File000507.nc File000508.nc File000509.nc File000510.nc File000511.nc File000512.nc File000513.nc File000514.nc File000515.nc File000516.nc File000517.nc File000518.nc File000519.nc File000520.nc File000521.nc File000522.nc File000523.nc File000524.nc File000525.nc File000526.nc File000527.nc File000528.nc File000529.nc File000530.nc File000531.nc File000532.nc File000533.nc File000534.nc File000535.nc File000536.nc File000537.nc File000538.nc File000539.nc File000540.nc File000541.nc File000542.nc File000543.nc File000544.nc File000545.nc File000546.nc File000547.nc File000548.nc File000549.nc File000550.nc File000551.nc File000552.nc File000553.nc File000554.nc File000555.nc File000556.nc File000557.nc File000558.nc File000559.nc File000560.nc File000561.nc File000562.nc File000563.nc File000564.nc File000565.nc File000566.nc File000567.nc File000568.nc File000569.nc File000570.nc File000571.nc File000572.nc File000573.nc File000574.nc File000575.nc IceVolumeHudson.nc
Tue Sep 13 22:41:13 2022: cdo -fldsum -mulc,10000 -mulc,10000 -selindexbox,279,477,230,362 -selvar,thk /scratch/m/m300792/HEDownload/HE66/HE66//pism_-064900/pism_-064900.nc Tmp/File000001.nc   CDO       @Climate Data Operators version 2.0.3 (https://mpimet.mpg.de/cdo)         time                standard_name         time   	long_name         time   units         seconds since 1-1-1    calendar      365_day    axis      T               -�   lon                standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X               -�   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y               -�   thk                       standard_name         land_ice_thickness     	long_name         land ice thickness     units         m      pism_intent       model_state             -�                �}Ȁ@�  C4�(w-o�}���   C4�yf(g��}�=   C4�0�>A��}�A�@  C4��r�n�}��9`  C3������}�·�  C3�L��|�}�5�  C2y &���}vC��  C1׻��T�}j�1�  C1aJ�w!C�}^İ   C1¦ݡ��}S.   C0�9���}GE�@  C0��	��};�*`  C0�dR�o�}/ƨ�  C0q�nG�}$&�  C0q�;�4�}G��  C0��l�k�}�"�  C0Ģ=ꍵ�} ȡ   C0��V{�|�	   C1�f���|�I�@  C13Q��#�|݊`  C1V���)�|�ʙ�  C1w�.{�e�|��  C1�+�uS�|�K��  C1���e\�|���  C1����
t�|�̒   C1�_G�Je�|�   C1׹�2��|�M�@  C1�W�6���|�`  C1��l�R=�|sΊ�  C2 ѓ#G�|h�  C2H쫹u�|\O��  C2:�p�s�|P��  C2]��?0�|DЃ   C2���/`�|9   C2�m>����|-Q@  C2�Cq�*�|!��`  C2�"<�PS�|�{�  C2�F�}Wx�|
��  C3�֠J�{�Sw�  C3�>�X��{���  C3!�7�{��t   C3,���o�{��   C37}����{�Up@  C3@�����{Õ�`  C3K�g�_�{��l�  C3a0<�xZ�{��  C3�C�����{�Wh�  C3�z� o��{����  C3�[�Ԇ�{��e   C3۳�v�{}�   C3�bMW���{qYa@  C4;�a2O�{e��`  C4�9d�*�{Y�]�  C4,�����{N۠  C4;D���"�{B[Y�  C4F��y�{6���  C4O#�fps�{*�V   C4V��ݏ�{�   C4]f�IO��{]R@  C4e�.{3S�{��`  C4wl9t�V�z��N�  C4���ڀ��z�̠  C4�H��h��z�_J�  C4�H:���z؟��  C4��w�b2�z��G   C4��a��E�z� �   C3[�#�N�z�aC@  C2��{3�z���`  C1�_�ڴP�z��?�  C1d������z�"��  C1Ej�F�z�c;�  C0Ъ�BP(�zz���  C0�\�i�]�zn�8   C0��*' �zc$�   C0j��N��zWe4@  C0[��2|6�zK��`  C0[�+^`��z?�0�  C0^�C��z4&��  C0s��]��z(g,�  C0��:2ݑ�z���  C0���8���z�)   C0��N���z(�   C1 ��nS�y�i%@  C1��`2p�y���`  C12!v��	�y��!�  C1Iq�߾�y�*��  C1V���F��y�k�  C1a�
�0�y����  C1o3�~�|�y��   C1~8V���y�,�   C1�9�AK�y�m@  C1���)n�y���`  C1��o&���y���  C2��>�\�yx.��  C2"�:tO�ylo�  C2@��m�y`���  C2Z�B~V�yT�   C2r��� ��yI0�   C2��x_�y=q@  C2��Aa�1�y1��`  C2�#	�J*�y%��  C2�x��7�y2��  C2�wu"{�yr��  C2�|��	�y�}�  C2����P�x���   C2�����V�x�4z   C3������x�t�@  C34dY����xӵv`  C3T	]4��x���  C3p��B ��x�6r�  C3���U���x�v��  C3��7J��x��n�  C3�ۄ�,f�x���   C3����x�8k   C3�ƣCd��x�x�@  C3�=�j��xu�g`  C3�R|���xi��  C3��}�:��x^:c�  C3��;�l]�xRz��  C3�/��1�xF�_�  C4D���w�x:��   C4.js��x/<\   C4J)�p���x#|�@  C4c���)�x�X`  C4|��X���x�ր  C4��^p�>�x >T�  C4�1��'�w�~��  C4�n�1LS�w�P�  C4������w���   C4����w�@M   C2���UU�wŀ�@  C2<�hU�t�w��I`  C1�U�-H�w�ǀ  C1!�]���w�BE�  C0���0�w����  C0��v��a�w��A�  C0��<��w�   C0ns��bv�wsD>   C0\�v���wg��@  C0S�}9��w[�:`  C0HX����wP��  C0@�բ�/�wDF6�  C0>%����w8���  C0:#�3t4�w,�2�  C0K���Gv�w!�   C0j9����wH/   C0��C��w	��@  C0��{n:��v��+`  C0�Ґ����v�	��  C0��&�R�v�J'�  C0Ϝ�̣�vڊ��  C0�Mֵ3��v��#�  C1�'�y��v��   C1>u�f��v�L    C1aNTS_�v���@  C1~����v��`  C1���G�v���  C1��T�6?�v�N�  C1ȯ��M8�v|���  C1��f��vp��  C1��2	5��ve�   C1�2ӈ�}�vYP   C2����9�vM��@  C2=g!�vA�`  C2 h�s���v6��  C27Հ�-��v*R	�  C2Y"Y.*�v���  C2�K�w���v��  C2�n����v�   C2�X�����u�T   C2���.@��u@  C2�χK��u���`  C3q�)��u�|�  C3 ����u�U��  C31�Xz��u��x�  C3?]��u��u����  C3J����u�u   C3Tn+~2��u�W�   C3]��Dh�u��q@  C3hh�,���u���`  C3|k~av?�uzm�  C3�^}��unY�  C3�/�u��ub�i�  C3ׇ�d�W�uV���  C3�>��P�uKf   C4��C��u?[�   C4!���u3�b@  C44�x69�u'��`  C4D@I2=��u^�  C4Q�iP�u]ܠ  C4\@��N�u�Z�  C4d"���t����  C4j^����t�W   C4o��G:��t�_�   C4w��,��tՠS@  C4�?�t�9�t���`  C4��8���t�!O�  C4���b���t�a͠  C4!�.�O��t��K�  C2��l�	��t����  C2X
��t�#H   C1vR?�t�c�   C0�b%Q���tw�D@  C0����tk��`  C0`�)~Ci�t`%@�  C00)U�/:�tTe��  C0�y|Za�tH�<�  C/�I�L> �t<��  C/��b/M�t1'9   C/m����t%g�   C/k]�`�n�t�5@  C/x��n�b�t�`  C/�b6yV�t)1�  C/��5<���s�i��  C/����	��s�-�  C/�u�3��s���  C/Έ��Af�s�+*   C0	ͽ���s�k�   C0'������s��&@  C0;���q��s��`  C0PA�����s�-"�  C0bGq�e
�s�m��  C0uL�걧�s���  C0��Fy�-�s���  C0����le�su/   C0����y]�sio�   C0�����s]�@  C0���ɯ_�sQ�`  C1�
h���sF1�  C1:�b���s:q��  C1Y!�����s.��  C1u)v�h�s"��  C1��"2���s3   C1��A;`��ss�   C1��w�r��@  C1��(D���r��`  C1�?F���r�5�  C1�M(��r�u��  C1��t��rж �  C2�5���r��~�  C2�x�h^�r�6�   C2<��8�r�w{   C2_ņ��O�r���@  C2�J����r��w`  C2�ш��r�8��  C2�͘W��r~ys�  C2�]��L��rr���  C2�'ұ�[�rf�o�  C3bpJf�r[:�   C3�j����rO{l   C3z��B��rC��@  C3*e����r7�h`  C35�'���r,<�  C3?/�-��r }d�  C3JT�r��r���  C3^�P�$_�r�`�  C3��?���q�>�   C3��'I���q�]   C3R홲�q��@  C3�LQ%�#�q� Y`  C3�F,�M��q�@׀  C4
ot���qU�  C4������q����  C4,�N ��q�Q�  C49Zl��q�B�   C4C���7��q��N   C4L<"[���q���@  C4R��c�q|J`  C4X�1e��qpDȀ  C4`���p�qd�F�  C4sJ4��s�qX���  C4��匷��qMB�  C4�m���qAF�   C4� �~���q5�?   C3dM4Y_Z�q)ǽ@  C2h�����q;`  C1�e�	���qH��  C1�t���q�7�  C0�L�j��p�ɵ�  C0j��L��p�
3�  C05 ���V�p�J�   C0
}R�ٷ�p׋0   C/�|����p�ˮ@  C/��N�n�p�,`  C/~���k��p�L��  C/f 9��Y�p��(�  C/k��G�p�ͦ�  C/`��L0M�p�$�  C/SVlr1��p�N�   C/Fq!�js�py�!   C/Bڜvd��pmϟ@  C/c�ŝ���pb`  C/�9�$��pVP��  C/��0�N��pJ��  C0W�'t�p>ї�  C0+����;�p3�  C0>�肝�p'R�   C0PS}�l_�p�   C0b�xi�S�pӐ@  C0u�p$�	�p`  C0�J}�z�o�   C0�� uh~�o�*@  C0�!\���o���  C0�'��Yy�o�,�  C1�Q�}�o��
   C1:�	����o{.@  C1Z���:X�oc��  C1t�&�1�oL/��  C1�w@"��o4��   C1�v�nD�o1�@  C1��ahs��o��  C1�ݟ"�n�3��  C1�������nִ�   C1���[c��n�5�@  C1�8����n���  C1�զ���n�7��  C2��B�nx��   C2,0����na9�@  C2M0�m�nI�Հ  C2j90�Q�n2;��  C2�ӄ��n��   C2���-&�n=�@  C2���ܕ/�m�ƀ  C2�ԩu��m�?��  C2�\���_�m���   C2⸩U��m�A�@  C2�V�r8�m�·�  C2�jfWn�mvC��  C3�V�5w�m^İ   C3~�l���mGE�@  C3%n�n��m/ƨ�  C3D�it���mG��  C3dxD.��m ȡ   C3�i�/S*�l�I�@  C3�7�m�i�l�ʙ�  C3�[E�L�l�K��  C3˟�M=~�l�̒   C3ު�>Q��l�M�@  C3��M����lsΊ�  C3���v� �l\O��  C4��@�H�lDЃ   C4�T�"�l-Q@  C4� �6-�l�{�  C4���D�k�Sw�  C4%�w�<<�k��t   C46pTO��k�Up@  C4LÈ�:<�k��l�  C3��GDQ��k�Wh�  C2���ӽ
�k��e   C1�y�	X9�kqYa@  C1-B9��kY�]�  C0��<0b�kB[Y�  C0d��sc�k*�V   C0*�ެ��k]R@  C/��ϥ ��j��N�  C/�y��p(�j�_J�  C/n��,���j��G   C/4�B+���j�aC@  C/�T
���j��?�  C.�#V,�j�c;�  C.��$�jn�8   C.���-m�jWe4@  C.����[��j?�0�  C.���V��j(g,�  C/����j�)   C/D�ʳ��i�i%@  C/���cd��i��!�  C/�%���s�i�k�  C0�xc��i��   C0^�%��i�m@  C05y�{�1�i���  C0E�=���ilo�  C0P݀����iT�   C0`D�ܓ3�i=q@  C0p�� �i%��  C0�Vl\���ir��  C0�!�s��h���   C0ԽĎX��h�t�@  C0������h���  C1F{!�g�h�v��  C19&�6���h���   C1UXW�17�h�x�@  C1n�F����hi��  C1��_[؋�hRz��  C1�ԗA���h:��   C1�K��?	�h#|�@  C1�����w�h�ր  C1��(���g�~��  C1ԧ�B"��g���   C1��"�$��gŀ�@  C1��'3��g�ǀ  C2�/���g����  C2?�W���g�   C2aER�i�gg��@  C2�JД��gP��  C2�
 �Q8�g8���  C2�FR8���g!�   C2�ё���g	��@  C2߲�g��f�	��  C2��Y����fڊ��  C2�Wf����f��   C3����f���@  C35�����f���  C3�<�M�f|���  C3*ѭr��fe�   C3?�u�fM��@  C3^Eg [�f6��  C3|�dp�f���  C3���/���f�   C3�#I�}��e@  C3�e�{;��e�|�  C3��ei��e��x�  C3�g[u�e�u   C3���(��e��q@  C32���d��ezm�  C2eQ#���eb�i�  C1��4���eKf   C1Ωh9f�e3�b@  C0� (�J��e^�  C0X��;��e�Z�  C0(Y�'���d�W   C0�v�d�dՠS@  C/��TK��d�!O�  C/؃=AN/�d��K�  C/�4kԍ2�d�#H   C/��G�4��dw�D@  C/�>WTa�d`%@�  C/�����\�dH�<�  C/�]��;�d1'9   C/�F߁�S�d�5@  C/�{>2���d)1�  C/��	�U��c�-�  C0ځ�c��c�+*   C0#n0���c��&@  C05{�����c�-"�  C0Q�D�V�c���  C0v��_��cu/   C0�$L:���c]�@  C0��El�L�cF1�  C0��d�i �c.��  C1J��]��c3   C1 ,�~��b��@  C19Q�z�a�b�5�  C1O.����bж �  C1b��]�b�6�   C1s�~���b���@  C1��t!�j�b�8��  C1�_�N#�br���  C1��ٔy��b[:�   C1��ظ��bC��@  C1�j="�b,<�  C1�tKh-M�b���  C2z���a�>�   C2'��s�e�a��@  C2F�ү���a�@׀  C2b��6�S�a����  C2{�i���a�B�   C2��OLVv�a���@  C2�'�y�9�apDȀ  C2�e�hG�aX���  C2Ź�����aAF�   C2�� �|G�a)ǽ@  C2�5����aH��  C2�R$D=��`�ɵ�  C2�)�o��`�J�   C3��r�`�ˮ@  C3'"L��`�`�L��  C3HՀFW�`�ͦ�  C3f�
9���`�N�   C3�C�\�`mϟ@  C3�9��',�`VP��  C3��P��4�`>ї�  C3���L��`'R�   C3�̓��^�`Ӑ@  C3�y4J��_�   C3�.D��_���  C3�`���_��
   C4�=i��_c��  C4#T����_4��   C4i���_��  C4&�ַE�^ִ�   C4B�ʋa�^���  C4^�]�� �^x��   C4tewS~�^I�Հ  C3�^��Yc�^��   C2�f'GVc�]�ƀ  C1�v��te�]���   C1��k@��]�·�  C0���v�!�]^İ   C0E�b�K4�]/ƨ�  C0
�ٯ8x�] ȡ   C/�v�=�u�\�ʙ�  C/��'��\�̒   C/P5��s�\sΊ�  C/)����,�\DЃ   C/%+��I�\�{�  C/4�G�#�[��t   C/I�HC���[��l�  C/h��IR��[��e   C/��6�[Y�]�  C/���1q�[*�V   C0�e����Z��N�  C09�w��Z��G   C0S߾,��Z��?�  C0i�3���Zn�8   C0|���-��Z?�0�  C0�� ���Z�)   C0��~�@A�Y��!�  C0��&��4�Y��   C0�B��B��Y���  C0ۖ�w��YT�   C0�s�7-�Y%��  C1�M7�X���   C1@F?'(��X���  C1`Bi��X���   C1}��]��Xi��  C1�ݯq�3�X:��   C1�_�'�c�X�ր  C1�eW�W���   C1ڒi���W�ǀ  C1�A�s<c�W�   C1���@<��WP��  C2�T�~��W!�   C2c���N�V�	��  C2��TD�V��   C27G;AN��V���  C2X�Y���Ve�   C2{
�v�"�V6��  C2���[�C�V�   C2����a�U�|�  C2���U�u   C2��{���Uzm�  C3�̉��UKf   C3LC�@�U^�  C3,m�WԶ�T�W   C3?�Ac��T�!O�  C3I��0��T�#H   C3O)�$#�T`%@�  C3V��%r��T1'9   C3ao-��T)1�  C3u��ҙ!�S�+*   C3�8p�i��S�-"�  C3�(vo�0�Su/   C3�Z�OP�SF1�  C3���n��S3   C4w���R�5�  C4������R�6�   C4"h묐=�R�8��  C3k.����R[:�   C2�Ĭ�s��R,<�  C1�a�"/"�Q�>�   C1,���p~�Q�@׀  C0���/�Q�B�   C0Jn�H�m�QpDȀ  C0K�;VG�QAF�   C/�>H���QH��  C/oz�v�P�J�   C/?2c��P�L��  C/�s_wX�P�N�   C.��{/���PVP��  C.����JA�P'R�   C.�=ǵ���O�   C.��R'�@�O��
   C.¾�J�d�O4��   C.�e&�/��Nִ�   C.�����Nx��   C.}�4���N��   C.���q�M���   C.��r�w�M^İ   C/k���(�M ȡ   C/=�fc�f�L�̒   C/�D#����LDЃ   C/����K��t   C0)�ee��K��e   C08M�a