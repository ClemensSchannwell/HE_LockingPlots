CDF  ?   
      time       lon       lat          	   CDI       @Climate Data Interface version 2.0.3 (https://mpimet.mpg.de/cdi)   Conventions       CF-1.5     source        #PISM development no-version-control    command      / /home/m/m300792/ModelCodes/mPISM_build/install_03-10-22/bin/pismr -i pism_-065000/pism_-065000.nc -stress_balance ssa+sia -y 100 -ys -65000 -o pism_-064900/pism_-064900.nc -surface_given_file forcing/smb_-064900.nc -surface given -ocean_th_file forcing/oce_-064900.nc -ocean th -topg_file pism_-065000/pism_-065000_new_topg.nc -bed_def prescribed -uplift_file pism_-065000/pism_-065000_dbdt.nc -extra_times 10 -extra_vars thk,usurf,bheatflx,velsurf_mag,velbase_mag,tillwat,bmelt,velbase,flux_divergence,velsurf,velbar,wvelsurf,wvelbase,climatic_mass_balance_cumulative,potential_climatic_mass_balance_cumulative,tempbase,tempsurf,mask,temppabase,tempicethk_basal,topg,nonneg_flux_cumulative,grounded_basal_flux_cumulative,floating_basal_flux_cumulative,discharge_flux_cumulative,Href,taud_mag,lat,lon,calving_mask,sea_level -extra_file pism_-064900/pism_-064900_extra -extra_split -ts_times yearly -ts_file pism_-064900/pism_-064900_ts.nc -topg_to_phi 15.0,30.0,-300.0,400.0 -ssafd_ksp_rtol 1e-7 -enthalpy_temperate_conductivity_ratio 0.001 -pseudo_plastic -pseudo_plastic_q .25 -pseudo_plastic_uthreshold 70 -till_effective_fraction_overburden 0.04 -tauc_slippery_grounding_lines -sia_e 8 -ssa_e 1 -ssa_upwind_factor 0.7 -calving eigen_calving,thickness_calving -thickness_calving_threshold 100 -part_redist -part_grid -cfbc -kill_icebergs -eigen_calving_K 1e16 -verbose 2 -sliding_scale_factor_reduces_tauc 1 -th_heat_coefficient 4e-5 -th_salinity_coefficient 4e-8 -o_format netcdf3 -sediment_file till_mask.nc -no_divide_tillphi -periodicity none -ocean_th_period 1 -o_order yxz
    proj4         j+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs     history_of_appended_files        �Fri Jul  8 17:07:29 2022: Appended file outFinal2 had following "history" attribute:
Fri Jul  8 17:07:27 2022: ncap2 -O -s lon_bnds=float(lon_bnds) outFinal1 outFinal2
Fri Jul  8 17:07:27 2022: ncap2 -O -s lat_bnds=float(lat_bnds) outFinal outFinal1
Fri Jul 08 17:07:25 2022: cdo -add thkNew thkOld outFinal
Fri Jul 08 17:07:24 2022: cdo -mul -selvar,thk pism_-064900/pism_-064900.nc MaskLaurentide thkNew
     NCO       _netCDF Operators version 5.0.6 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)    history       |Thu Sep 15 09:18:02 2022: cdo -mergetime File000001.nc File000002.nc File000003.nc File000004.nc File000005.nc File000006.nc File000007.nc File000008.nc File000009.nc File000010.nc File000011.nc File000012.nc File000013.nc File000014.nc File000015.nc File000016.nc File000017.nc File000018.nc File000019.nc File000020.nc File000021.nc File000022.nc File000023.nc File000024.nc File000025.nc File000026.nc File000027.nc File000028.nc File000029.nc File000030.nc File000031.nc File000032.nc File000033.nc File000034.nc File000035.nc File000036.nc File000037.nc File000038.nc File000039.nc File000040.nc File000041.nc File000042.nc File000043.nc File000044.nc File000045.nc File000046.nc File000047.nc File000048.nc File000049.nc File000050.nc File000051.nc File000052.nc File000053.nc File000054.nc File000055.nc File000056.nc File000057.nc File000058.nc File000059.nc File000060.nc File000061.nc File000062.nc File000063.nc File000064.nc File000065.nc File000066.nc File000067.nc File000068.nc File000069.nc File000070.nc File000071.nc File000072.nc File000073.nc File000074.nc File000075.nc File000076.nc File000077.nc File000078.nc File000079.nc File000080.nc File000081.nc File000082.nc File000083.nc File000084.nc File000085.nc File000086.nc File000087.nc File000088.nc File000089.nc File000090.nc File000091.nc File000092.nc File000093.nc File000094.nc File000095.nc File000096.nc File000097.nc File000098.nc File000099.nc File000100.nc File000101.nc File000102.nc File000103.nc File000104.nc File000105.nc File000106.nc File000107.nc File000108.nc File000109.nc File000110.nc File000111.nc File000112.nc File000113.nc File000114.nc File000115.nc File000116.nc File000117.nc File000118.nc File000119.nc File000120.nc File000121.nc File000122.nc File000123.nc File000124.nc File000125.nc File000126.nc File000127.nc File000128.nc File000129.nc File000130.nc File000131.nc File000132.nc File000133.nc File000134.nc File000135.nc File000136.nc File000137.nc File000138.nc File000139.nc File000140.nc File000141.nc File000142.nc File000143.nc File000144.nc File000145.nc File000146.nc File000147.nc File000148.nc File000149.nc File000150.nc File000151.nc File000152.nc File000153.nc File000154.nc File000155.nc File000156.nc File000157.nc File000158.nc File000159.nc File000160.nc File000161.nc File000162.nc File000163.nc File000164.nc File000165.nc File000166.nc File000167.nc File000168.nc File000169.nc File000170.nc File000171.nc File000172.nc File000173.nc File000174.nc File000175.nc File000176.nc File000177.nc File000178.nc File000179.nc File000180.nc File000181.nc File000182.nc File000183.nc File000184.nc File000185.nc File000186.nc File000187.nc File000188.nc File000189.nc File000190.nc File000191.nc File000192.nc File000193.nc File000194.nc File000195.nc File000196.nc File000197.nc File000198.nc File000199.nc File000200.nc File000201.nc File000202.nc File000203.nc File000204.nc File000205.nc File000206.nc File000207.nc File000208.nc File000209.nc File000210.nc File000211.nc File000212.nc File000213.nc File000214.nc File000215.nc File000216.nc File000217.nc File000218.nc File000219.nc File000220.nc File000221.nc File000222.nc File000223.nc File000224.nc File000225.nc File000226.nc File000227.nc File000228.nc File000229.nc File000230.nc File000231.nc File000232.nc File000233.nc File000234.nc File000235.nc File000236.nc File000237.nc File000238.nc File000239.nc File000240.nc File000241.nc File000242.nc File000243.nc File000244.nc File000245.nc File000246.nc File000247.nc File000248.nc File000249.nc File000250.nc File000251.nc File000252.nc File000253.nc File000254.nc File000255.nc File000256.nc File000257.nc File000258.nc File000259.nc File000260.nc File000261.nc File000262.nc File000263.nc File000264.nc File000265.nc File000266.nc File000267.nc File000268.nc File000269.nc File000270.nc File000271.nc File000272.nc File000273.nc File000274.nc File000275.nc File000276.nc File000277.nc File000278.nc File000279.nc File000280.nc File000281.nc File000282.nc File000283.nc File000284.nc File000285.nc File000286.nc File000287.nc File000288.nc File000289.nc File000290.nc File000291.nc File000292.nc File000293.nc File000294.nc File000295.nc File000296.nc File000297.nc File000298.nc File000299.nc File000300.nc File000301.nc File000302.nc File000303.nc File000304.nc File000305.nc File000306.nc File000307.nc File000308.nc File000309.nc File000310.nc File000311.nc File000312.nc File000313.nc File000314.nc File000315.nc File000316.nc File000317.nc File000318.nc File000319.nc File000320.nc File000321.nc File000322.nc File000323.nc File000324.nc File000325.nc File000326.nc File000327.nc File000328.nc File000329.nc File000330.nc File000331.nc File000332.nc File000333.nc File000334.nc File000335.nc File000336.nc File000337.nc File000338.nc File000339.nc File000340.nc File000341.nc File000342.nc File000343.nc File000344.nc File000345.nc File000346.nc File000347.nc File000348.nc File000349.nc File000350.nc File000351.nc File000352.nc File000353.nc File000354.nc File000355.nc File000356.nc File000357.nc File000358.nc File000359.nc File000360.nc File000361.nc File000362.nc File000363.nc File000364.nc File000365.nc File000366.nc File000367.nc File000368.nc File000369.nc File000370.nc File000371.nc File000372.nc File000373.nc File000374.nc File000375.nc File000376.nc File000377.nc File000378.nc File000379.nc File000380.nc File000381.nc File000382.nc File000383.nc File000384.nc File000385.nc File000386.nc File000387.nc File000388.nc File000389.nc File000390.nc File000391.nc File000392.nc File000393.nc File000394.nc File000395.nc File000396.nc File000397.nc File000398.nc File000399.nc File000400.nc File000401.nc File000402.nc File000403.nc File000404.nc File000405.nc File000406.nc File000407.nc File000408.nc File000409.nc File000410.nc File000411.nc File000412.nc File000413.nc File000414.nc File000415.nc File000416.nc File000417.nc File000418.nc File000419.nc File000420.nc File000421.nc File000422.nc File000423.nc File000424.nc File000425.nc File000426.nc File000427.nc File000428.nc File000429.nc File000430.nc File000431.nc File000432.nc File000433.nc File000434.nc File000435.nc File000436.nc File000437.nc File000438.nc File000439.nc File000440.nc File000441.nc File000442.nc File000443.nc File000444.nc File000445.nc File000446.nc File000447.nc File000448.nc File000449.nc File000450.nc File000451.nc File000452.nc File000453.nc File000454.nc File000455.nc File000456.nc File000457.nc File000458.nc File000459.nc File000460.nc File000461.nc File000462.nc File000463.nc File000464.nc File000465.nc File000466.nc File000467.nc File000468.nc File000469.nc File000470.nc File000471.nc File000472.nc File000473.nc File000474.nc File000475.nc File000476.nc File000477.nc File000478.nc File000479.nc File000480.nc File000481.nc File000482.nc File000483.nc File000484.nc File000485.nc File000486.nc File000487.nc File000488.nc File000489.nc File000490.nc File000491.nc File000492.nc File000493.nc File000494.nc File000495.nc File000496.nc File000497.nc File000498.nc File000499.nc File000500.nc File000501.nc File000502.nc File000503.nc File000504.nc File000505.nc File000506.nc File000507.nc File000508.nc File000509.nc File000510.nc File000511.nc File000512.nc File000513.nc File000514.nc File000515.nc File000516.nc File000517.nc File000518.nc File000519.nc File000520.nc File000521.nc File000522.nc File000523.nc File000524.nc File000525.nc File000526.nc File000527.nc File000528.nc File000529.nc File000530.nc File000531.nc File000532.nc File000533.nc File000534.nc File000535.nc File000536.nc File000537.nc File000538.nc File000539.nc File000540.nc File000541.nc File000542.nc File000543.nc File000544.nc File000545.nc File000546.nc File000547.nc File000548.nc File000549.nc File000550.nc File000551.nc File000552.nc File000553.nc File000554.nc File000555.nc File000556.nc File000557.nc File000558.nc File000559.nc File000560.nc File000561.nc File000562.nc File000563.nc File000564.nc File000565.nc File000566.nc File000567.nc File000568.nc File000569.nc File000570.nc File000571.nc File000572.nc File000573.nc File000574.nc File000575.nc IceVolumeKenzie.nc
Thu Sep 15 09:13:17 2022: cdo -fldsum -mulc,10000 -mulc,10000 -selindexbox,211,427,600,415 -selvar,thk /work/ba0989/m300792/HE_Runs/ExperimentsComposite//HE75//pism_-064900/pism_-064900.nc Tmp/File000001.nc   CDO       @Climate Data Operators version 2.0.3 (https://mpimet.mpg.de/cdo)         time                standard_name         time   	long_name         time   units         seconds since 1-1-1    calendar      365_day    axis      T               -�   lon                standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X               -�   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y               -�   thk                       standard_name         land_ice_thickness     	long_name         land ice thickness     units         m      pism_intent       model_state             -�                �}Ȁ@�  C-T�I�T�}���   C-�G[KG��}�=   C-��{Ԃ��}�A�@  C-�?i}��}��9`  C.�Դ�7�}�·�  C.Ih����}�5�  C.y����%�}vC��  C.�mxn5�}j�1�  C.�<Re�=�}^İ   C.�\����}S.   C/�eR���}GE�@  C/'�6k�>�};�*`  C/@x{�/�}/ƨ�  C/]��$+5�}$&�  C/{����}G��  C/��X(l(�}�"�  C/��m
���} ȡ   C/���U��|�	   C/�_�j���|�I�@  C/�x��\��|݊`  C0zLalV�|�ʙ�  C0].Cv�k�|��  C0�*�X1�|�K��  C0�H�3��|���  C0���"���|�̒   C0�=?�c��|�   C1/���|�M�@  C1� uW�|�`  C1,�1P@O�|sΊ�  C1>�;�;`�|h�  C1K�h�Qp�|\O��  C1R1����|P��  C1Y�����|DЃ   C1_��<���|9   C1c�,� ��|-Q@  C1g�Hx��|!��`  C0�qb��|�{�  C/�M��**�|
��  C.��F���{�Sw�  C.p�U�o�{���  C-�5�ϋ�{��t   C-]�G�/��{��   C-��ܸn�{�Up@  C,�`�c�{Õ�`  C,�줜��{��l�  C,w"�����{��  C,S��V���{�Wh�  C,.������{����  C,$�xmq��{��e   C,"nK��{}�   C,!� ����{qYa@  C,*��\��{e��`  C,I��x&�{Y�]�  C,cYdaO�{N۠  C,u�S�I��{B[Y�  C,��j�{6���  C,��T5ۏ�{*�V   C,�$تY��{�   C,��HE�{]R@  C,����s�{��`  C-\�2g�z��N�  C-q��ǆ1�z�̠  C-��K��z�_J�  C-���i�w�z؟��  C-�[����z��G   C._��r�z� �   C.?����z�aC@  C.o�����z���`  C.��sײ	�z��?�  C.9Kv��z�"��  C.�<P6�>�z�c;�  C/�7=o�zz���  C/c���8�zn�8   C/�m��L��zc$�   C/�Y� 7�zWe4@  C0b����zK��`  C0\/����z?�0�  C0(
�O6�z4&��  C01u^�\o�z(g,�  C0@��C�+�z���  C0I؜Va3�z�)   C0`DKa ��z(�   C0����t�y�i%@  C0�����q�y���`  C0��3H-��y��!�  C0�m�m�y�*��  C0�����U�y�k�  C0�cń���y����  C1������y��   C1�`xl��y�,�   C/�=ڜR��y�m@  C.���q��y���`  C-ҨQ����y���  C-J6�����yx.��  C,�Ts����ylo�  C,�|��GU�y`���  C,a2�$��yT�   C,T����yI0�   C,aѭ�|F�y=q@  C,ne^��p�y1��`  C,���H�y%��  C,��<�X�y2��  C,�qپ���yr��  C-	T�jU�y�}�  C-1���}��x���   C-i\��5��x�4z   C-����#��x�t�@  C-���3��xӵv`  C.��jA��x���  C.;_t�i�x�6r�  C.v먨���x�v��  C.�ʏa���x��n�  C.�hT�*��x���   C/�Y�$�x�8k   C/��9v��x�x�@  C/9�+��f�xu�g`  C/P���޹�xi��  C/ciR��x^:c�  C/m����xRz��  C/P	�6���xF�_�  C/[d�">�x:��   C/��>T�'�x/<\   C/�2LF��x#|�@  C0�r�<��x�X`  C07���w�x�ր  C0W #U�m�x >T�  C0uB�RO��w�~��  C0�#�H��w�P�  C0���Cf�w���   C0�.�H!�w�@M   C0׽��S�wŀ�@  C0��I��w��I`  C0�~5�CK�w�ǀ  C0�n}C�e�w�BE�  C1�c��w����  C1�`���w��A�  C1G�1k��w�   C10L�_$�wsD>   C1p1 �wg��@  C1"���]�w[�:`  C1A;�����wP��  C1\�re���wDF6�  C1k�Sx}�w8���  C1�oiy-m�w,�2�  C1�[{V� �w!�   C1��
��g�wH/   C11�6dx�w	��@  C/P����U�v��+`  C-�3���v�	��  C-������v�J'�  C,x���a�vڊ��  C,+�����v��#�  C,+��(��v��   C,T=u���v�L    C,c���S�v���@  C,n]մ_�v��`  C,u��Ֆ,�v���  C,|01"���v�N�  C,�K�;���v|���  C,�L�nR�vp��  C,ཾ�QK�ve�   C-G<��V�vYP   C-�N�{�,�vM��@  C."���d�vA�`  C.w��]��v6��  C-���.:��v*R	�  C.�v��v���  C.�Q��v��  C.?�nR��v�   C.dZ���E�u�T   C.��;3���u@  C.��엛��u���`  C.�ur~Pq�u�|�  C.��D��u�U��  C/�.9�u��x�  C/+�6��u����  C/Q�n�j�u�u   C/x(����u�W�   C/����`;�u��q@  C/����u���`  C0&�g���uzm�  C0G��u��unY�  C0V���=�ub�i�  C0p�?���uV���  C0�ŨDM�uKf   C0��0ʧ��u?[�   C0�c,C��u3�b@  C0�h� ��u'��`  C0�u	�h��u^�  C19�8�u]ܠ  C1 ���!��u�Z�  C1'	"���t����  C0�N�<��t�W   C/��8���t�_�   C.�zX1�k�tՠS@  C.
j��`��t���`  C-r�6$�l�t�!O�  C,��,��J�t�a͠  C,��v�؃�t��K�  C,q��=�A�t����  C,����;�t�#H   C,�6t7A��t�c�   C,�7Ϸ��tw�D@  C,�e����tk��`  C,�2o��.�t`%@�  C,����m�tTe��  C,�0����tH�<�  C,����_�t<��  C,��zp��t1'9   C-Q���t%g�   C-3��g�t�5@  C-Y�H�a�t�`  C-{��~@��t)1�  C-�K�����s�i��  C-�IKw�|�s�-�  C-���ŷ�s���  C-�����s�+*   C-��X~!�s�k�   C.
'm�v�s��&@  C.*�X�
t�s��`  C.r�=AD�s�-"�  C.�3Fǀ_�s�m��  C.�?霡�s���  C/��d��s���  C/M�Z��su/   C/�zi]�Y�sio�   C/���[�y�s]�@  C0��5�}�sQ�`  C0"#�R���sF1�  C06ս�kP�s:q��  C/�~]���s.��  C.�ծB�c�s"��  C.0�륋��s3   C-���^�ss�   C-{*�!��r��@  C-2�I��Z�r��`  C,���F���r�5�  C,����cT�r�u��  C,��}����rж �  C,z�.H�1�r��~�  C,��d"�r�6�   C,��ʵ@q�r�w{   C,�֌����r���@  C,�5`��r��w`  C,����Pi�r�8��  C,�L�2���r~ys�  C,NI��!��rr���  C,9Ǫ.8�rf�o�  C,.�99�r[:�   C,1��w��rO{l   C,@�U����rC��@  C,CS���u�r7�h`  C,O�����r,<�  C,f\*���r }d�  C,��%6!�r���  C,�Z�GNf�r�`�  C,��P�|�q�>�   C,��J"���q�]   C,��J� V�q��@  C,��w�DE�q� Y`  C-$'�B��q�@׀  C-l�$k���qU�  C-r��=���q����  C-�ed���q�Q�  C-�W���q�B�   C-�(H�b��q��N   C.3w��JA�q���@  C.f(�w�q|J`  C.��p��qpDȀ  C.�7QD{P�qd�F�  C.�Jҁ��qX���  C/
x*@�x�qMB�  C/$�ۛ�qAF�   C/9�� ��q5�?   C/T�Z����q)ǽ@  C/kpIM~�q;`  C/4J�� �qH��  C/��)6�q�7�  C/ΗԈ�x�p�ɵ�  C0`RpZw�p�
3�  C07���h�p�J�   C0w&\l���p׋0   C0�D|����p�ˮ@  C0��0T�p�,`  C0��x���p�L��  C0ʤ�����p��(�  C/<�S���p�ͦ�  C-�}7�q�p�$�  C-�@I�N�p�N�   C,��o4��py�!   C,MNzO�pmϟ@  C,
I�d=�pb`  C+�uAT�pVP��  C+��*���pJ��  C+�P��U�p>ї�  C+�����p3�  C+�˜��p'R�   C+�B)il�p�   C+���@���pӐ@  C+�9p)���p`  C+��&�c��o�   C,W�~���o�*@  C,�Z!'g��o���  C-:u�I*>�o�,�  C-��k!r��o��
   C-�q���c�o{.@  C-�����,�oc��  C-�vN���oL/��  C-ݽL9և�o4��   C. �#�c��o1�@  C.'�ԛ���o��  C.K_Q|:�n�3��  C.fa~}j��nִ�   C.}�B�F��n�5�@  C.��VK]�n���  C.�"�QM6�n�7��  C.����nx��   C.핍Ĉ
�na9�@  C/!��^�nI�Հ  C/k�ljg�n2;��  C/������n��   C0+�ҍ���n=�@  C0n�L�Y�m�ƀ  C0z^���A�m�?��  C0����A�m���   C0�o�q���m�A�@  C0������m�·�  C0����Q�mvC��  C0�oYXW�m^İ   C0���Z7�mGE�@  C0�1�X��m/ƨ�  C/jr��}��mG��  C.��W2_�m ȡ   C-����]�l�I�@  C-�����l�ʙ�  C,�ӵ,��l�K��  C,�g�}��l�̒   C+���BN��l�M�@  C+�?_��j�lsΊ�  C+|zC�9�l\O��  C+��5���lDЃ   C,?�O�l-Q@  C,��:1N��l�{�  C,ߒ��H�k�Sw�  C-;饡���k��t   C-�e��k��k�Up@  C-�~s�%�k��l�  C-�{�����k�Wh�  C-�Oj��8�k��e   C-�N�Ó��kqYa@  C-��n��kY�]�  C.m�Z���kB[Y�  C.9��DP�k*�V   C.]����g�k]R@  C.}�q�+b�j��N�  C.�� �o��j�_J�  C.��oP��j��G   C.��" �r�j�aC@  C.ѿA�#�j��?�  C.��D?��j�c;�  C/1�¾��jn�8   C/��m۹)�jWe4@  C/�S��^�j?�0�  C/�O 	j�j(g,�  C0���c�j�)   C0(�X�e��i�i%@  C0B%��7;�i��!�  C0Y����i�k�  C0pc�G:�i��   C0��X����i�m@  C0��J )�i���  C0��8xff�ilo�  C0�Jp~���iT�   C0�Ɋ�i�i=q@  C1�{�ݻ�i%��  C1�7���ir��  C1!��;�h���   C1 J;8�h�t�@  C0����2�h���  C0�<�I[��h�v��  C/`�9�t�h���   C.OM���(�h�x�@  C-��,����hi��  C-�;�[�hRz��  C,�:���0�h:��   C,u��J��h#|�@  C,AzJ{ �h�ր  C,$�kK{��g�~��  C,O�C��g���   C,:2WS�gŀ�@  C, x�k�g�ǀ  C,�u1��g����  C,V1�{�g�   C,8+�(j�gg��@  C,U&��
��gP��  C,h1�*���g8���  C,p�&N��g!�   C,y?8����g	��@  C,�זց��f�	��  C,�L/���fڊ��  C,���UH�f��   C-G��h)��f���@  C-�\z�6��f���  C.I�
[���f|���  C.o�Q�A�fe�   C.L�[h��fM��@  C.J��L>�f6��  C.D�p3���f���  C.T�T�)p�f�   C.iu_ڸ��e@  C.���h���e�|�  C.����!�e��x�  C.ٕ�aRQ�e�u   C/P�]D�e��q@  C/\EC���ezm�  C/�O��S2�eb�i�  C/����|��eKf   C/�.#����e3�b@  C/���PH��e^�  C0
fK�	�e�Z�  C0)
\��.�d�W   C06:�/���dՠS@  C0ROȿ�V�d�!O�  C0w4s����d��K�  C0�d*5_��d�#H   C0�UA� �dw�D@  C0
�U�Tv�d`%@�  C.ā~��K�dH�<�  C.#�Ԏ�d1'9   C-�6���2�d�5@  C-0���d)1�  C,�K�%#L�c�-�  C,��#��]�c�+*   C,U\��e��c��&@  C,?�n��#�c�-"�  C,5��c�+�c���  C,0���cu/   C,?d'�ѱ�c]�@  C,H�H�a�cF1�  C,[�e�BJ�c.��  C,�wĽC��c3   C,�Q��I[�b��@  C-\e�e���b�5�  C-L�}-��bж �  C-l��A.��b�6�   C-���$��b���@  C-�'�0�b�8��  C-�<�]��br���  C-�3c��b[:�   C.�S�Z�bC��@  C.(l��D�b,<�  C.W
�4�P�b���  C.���{0o�a�>�   C.��. *�a��@  C.�[����a�@׀  C.���1H�a����  C.��E���a�B�   C/�A�7�a���@  C/�zqP��apDȀ  C/00&(�aX���  C/I|����aAF�   C/�6
�9�a)ǽ@  C0�Ls���aH��  C02�8vN��`�ɵ�  C0c�S&�`�J�   C0��d��%�`�ˮ@  C0����&_�`�L��  C0���0��`�ͦ�  C0�b��p��`�N�   C0�eEuP��`mϟ@  C0��sTJ��`VP��  C0�wz�mP�`>ї�  C1R��=��`'R�   C000�,���`Ӑ@  C.��k4+
�_�   C.;�����_���  C-�
�O>��_��
   C,�28͗C�_c��  C,��h���_4��   C,Y�!����_��  C,ha4�O<�^ִ�   C,�4O� ��^���  C--�M�2�^x��   C-�f]��u�^I�Հ  C-�e�ޝJ�^��   C-��9'���]�ƀ  C-�a���-�]���   C-������]�·�  C-���-~-�]^İ   C-�~ya�(�]/ƨ�  C-�"ָƔ�] ȡ   C.�h崞�\�ʙ�  C.,�q�'6�\�̒   C.Q�����\sΊ�  C.q2d/+z�\DЃ   C.�VWzJ'�\�{�  C.���oP)�[��t   C.�����i�[��l�  C.������[��e   C.�QJuJ�[Y�]�  C/�YS�n�[*�V   C/9���K�Z��N�  C/[�4Y��Z��G   C/��OK���Z��?�  C/��<�2��Zn�8   C0) Y���Z?�0�  C0W(���Z�)   C0y�:���Y��!�  C0��p	E<�Y��   C0���µ��Y���  C0ĳ.3�u�YT�   C0��uQ�|�Y%��  C0�ұ��X���   C0�3����X���  C.�r?M&�X���   C-��Y7���Xi��  C-d�nf��X:��   C,ܭm%�4�X�ր  C,r��]��W���   C,&p��wu�W�ǀ  C,���W�   C+�P���B�WP��  C,	���R�W!�   C,L�0�V�	��  C,&����I�V��   C,*G�q�:�V���  C,��!k��Ve�   C+��&�a�V6��  C+�����\�V�   C+�"	V��U�|�  C+���b���U�u   C+�8|�ʅ�Uzm�  C+zc9���UKf   C+k�r����U^�  C+_��'��T�W   C+j�6�G��T�!O�  C+�1`M�T�#H   C+�b�
�T`%@�  C+��и�T1'9   C+�:q��T)1�  C,�?Ė5�S�+*   C,���#�S�-"�  C-�(�O�Su/   C-tGP�z�SF1�  C-�� �P�S3   C-��v)5��R�5�  C.¶4���R�6�   C.=�JFOI�R�8��  C.u�_	�$�R[:�   C.�̪h�u�R,<�  C.�f��u��Q�>�   C.��?f���Q�@׀  C/�S���Q�B�   C/E�n��y�QpDȀ  C/i�̺%��QAF�   C/�H���QH��  C/��#��d�P�J�   C/��(�ˋ�P�L��  C/Ҙ-�\�P�N�   C/�����PVP��  C0 o)D}J�P'R�   C0T����k�O�   C0d �g�O��
   C0}�f����O4��   C0�[�pL��Nִ�   C00���q��Nx��   C.�<<�=�N��   C-����4��M���   C-)����M^İ   C,�R��.R�M ȡ   C,B�x�؅�L�̒   C+�X|V�)�LDЃ   C+�ۅ)���K��t   C+r��
=��K��e   C+C�&j-k