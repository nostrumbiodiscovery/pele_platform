import pytest
import sys
import os
import subprocess
from test_config import parse_template
sys.path.insert(1, os.path.join(sys.path[0], '..'))
import PlpRotTemp.PlopRotTemp as pl
import filecmp
#from ciphers.DESchiper import convert_to_binary, DES, XOR


MAE_FILE = 'ain.mae'
ROOT = 'ain'
MAE_CONVERSION = [0,1,2,3,4,13,5,14,6,15,7,16,8,9,10,11,12,17,18,19]
TEMPLATE_CONVERSION = [0,1,2,3,4,6,8,10,12,13,14,15,16,5,7,9,11,17,18,19]
BONDS = [[0, 1], [1, 2], [1, 3], [3, 4], [3, 8], [4, 5], [4, 13], [5, 6], 
        [5, 14], [6, 7], [6, 15], [7, 8], [7, 16], [8, 9], [9, 10], [10, 11],
        [10, 12], [12, 17], [12, 18], [12, 19]]
PARENT_RESULT = [0, 1, 2, 2, 4, 5, 6, 7, 4, 9, 10, 11, 11, 5, 6, 7, 8, 13, 13, 13]
REPOS_PATH = os.path.abspath(os.path.join(__file__ ,"../.."))
TEST_PATH = os.path.join(os.path.dirname(__file__), 'data')
MAIN_PATH =  os.path.join(REPOS_PATH, 'PlpRotTemp/main.py')
OLD_MAIN_PATH = '/home/dani/repos/presentation/PlpRotTemp/main_16_9_17.py'
try:
    PYTHON_PATH = os.path.join(os.environ['SCHRODINGER'] + "/utilities/python")
except KeyError:
    print("Set SCHRODINGER environment variable path")



# @pytest.mark.parametrize("argument, expected", [
#                          ('ain.mae', 'ain'),
#                          ('~/ain.mae', 'ain'),
#                          ('/opt/schrodinger/ain.mae', 'ain'),
#                          ])
# def test_get_root(argument, expected):
#     root = pl.get_root_path(argument)
#     assert root == expected

# # @pytest.mark.parametrize("MAE_FILE, root, OPLS, hetgrp_opt, old_name, new_name, expected", [
# #                          (os.path.join(TEST_PATH, 'ain.mae'), 'ain', '2005', '', '', '', 'ain.hetgrp_ffgen'),
# #                          (os.path.join(TEST_PATH,'MI4.mae'), 'MI4', '2005', '', '', '', 'mi4.hetgrp_ffgen'),
# #                          ])
# # def test_build_template(MAE_FILE, root, OPLS, hetgrp_opt, old_name, new_name, expected):
# #     [template_file, output_template_file, mae_file_hetgrp_ffgen, files, resname] = pl.build_template(MAE_FILE, root, OPLS, hetgrp_opt, old_name, new_name)
# #     assert template_file == expected




# @pytest.mark.parametrize("mae_file, template_file, mae_expected, template_expected", [
#                          (os.path.join(TEST_PATH, 'ain.mae'), 'ain.hetgrp_ffgen' , MAE_CONVERSION, TEMPLATE_CONVERSION),
#                          ])
# def test_MatchTempMaeAtoms(mae_file, template_file, mae_expected, template_expected):
#     [mae2temp, temp2mae] = pl.MatchTempMaeAtoms(mae_file, template_file)
#     assert mae2temp == template_expected
#     assert temp2mae == mae_expected

# @pytest.mark.parametrize("mae_file, expected", [
#                          (os.path.join(TEST_PATH, 'ain.mae'), ''),
#                          (os.path.join(TEST_PATH, 'ain_repited.mae'), Exception),
#                          ])
# def test_check_repite_names(mae_file, expected):
#     atomnames = pl.find_names_in_mae(mae_file)
#     if(expected == Exception):
#         with pytest.raises(expected):
#             pl.check_repite_names(atomnames)
#     else:
#         pass



# @pytest.mark.parametrize("rotamer_library", [
#                          (os.path.join(TEST_PATH,'ain_vdw')),
#                          ])
# def test_check_replace_vdwr(rotamer_library):
#     pl.replace_vdwr_from_library(rotamer_library)
#     radius_vdw_info, start_index, end_index = pl.parse_nonbonded(rotamer_library)
#     for i, rdw_line in enumerate(radius_vdw_info):
#         NBOND_info = rdw_line.split()
#         rdw = float(NBOND_info[1])/2.0
#         if(rdw== 0):
#             assert 

 

# @pytest.mark.parametrize("input_file", [
                         #   (MAE_FILE),
                         #  ('MI4.mae'),
                         #  ('1o3p_ligand.mae'),
                         #  ('1.mae'),('2.mae'),('3.mae'),('4.mae'),('5.mae'),('6.mae'),('7.mae'),('8.mae'),('9.mae'),('10.mae'),
                         #  ('11.mae'),('12.mae'),('13.mae'),('14.mae'),('15.mae'),('16.mae'),('17.mae'),('18.mae'),('19.mae'),('20.mae'),
                         #  ('21.mae'),('22.mae'),('23.mae'),('24.mae'),('25.mae'),('26.mae'),('27.mae'),('28.mae'),('29.mae'),('30.mae'),
                         #  ('31.mae'),('32.mae'),('33.mae'),('34.mae'),('35.mae'),('36.mae'),('37.mae'),('38.mae'),('39.mae'),('40.mae'),
                         #  ('41.mae'),('42.mae'),('43.mae'),('44.mae'),('45.mae'),('46.mae'),('47.mae'),('48.mae'),('49.mae'),('50.mae'),
                         #  ('51.mae'),('52.mae'),('53.mae'),('54.mae'),('55.mae'),('56.mae'),('57.mae'),('58.mae'),('59.mae'),('60.mae'),
                         #  ('61.mae'),('62.mae'),('63.mae'),('64.mae'),('65.mae'),('66.mae'),('67.mae'),('68.mae'),('69.mae'),('70.mae'),
                         #  ('71.mae'),('72.mae'),('73.mae'),('74.mae'),('75.mae'),('76.mae'),('77.mae'),('78.mae'),('79.mae'),('80.mae'),
                         #  ('81.mae'),('82.mae'),('83.mae'),('84.mae'),('85.mae'),('86.mae'),('87.mae'),('88.mae'),('89.mae'),('90.mae'),
                         #  ('91.mae'),('92.mae'),('93.mae'),('94.mae'),('95.mae'),('96.mae'),('97.mae'),('98.mae'),('99.mae'),('100.mae'),

                         #  ('100.mae'),('101.mae'),('102.mae'),('103.mae'),('104.mae'),('105.mae'),('106.mae'),('107.mae'),('108.mae'),('109.mae'),('110.mae'),
                         #  ('110.mae'),('111.mae'),('112.mae'),('113.mae'),('114.mae'),('115.mae'),('116.mae'),('117.mae'),('118.mae'),('119.mae'),('120.mae'),
                         #  ('120.mae'),('121.mae'),('122.mae'),('123.mae'),('124.mae'),('125.mae'),('126.mae'),('127.mae'),('128.mae'),('129.mae'),('130.mae'),
                         #  ('131.mae'),('132.mae'),('133.mae'),('134.mae'),('135.mae'),('136.mae'),('137.mae'),('138.mae'),('139.mae'),('140.mae'),
                         #  ('141.mae'),('142.mae'),('143.mae'),('144.mae'),('145.mae'),('146.mae'),('147.mae'),('148.mae'),('149.mae'),('150.mae'),
                         #  ('151.mae'),('152.mae'),('153.mae'),('154.mae'),('155.mae'),('156.mae'),('157.mae'),('158.mae'),('159.mae'),('160.mae'),
                         #  ('161.mae'),('162.mae'),('163.mae'),('164.mae'),('165.mae'),('166.mae'),('167.mae'),('168.mae'),('169.mae'),('170.mae'),
                         #  ('171.mae'),('172.mae'),('173.mae'),('174.mae'),('175.mae'),('176.mae'),('177.mae'),('178.mae'),('179.mae'),('180.mae'),
                         #  ('181.mae'),('182.mae'),('183.mae'),('184.mae'),('185.mae'),('186.mae'),('187.mae'),('188.mae'),('189.mae'),('190.mae'),
                         #  ('191.mae'),('192.mae'),('193.mae'),('194.mae'),('195.mae'),('196.mae'),('197.mae'),('198.mae'),('199.mae'),('200.mae'),

                         #  ('201.mae'), ('202.mae'), ('203.mae'), ('204.mae'), ('205.mae'), ('206.mae'), ('207.mae'), ('208.mae'), ('209.mae'), ('210.mae'),
                         #  ('211.mae'), ('212.mae'), ('213.mae'), ('214.mae'), ('215.mae'), ('216.mae'), ('217.mae'), ('218.mae'), ('219.mae'), ('220.mae'),
                         #  ('221.mae'), ('222.mae'), ('223.mae'), ('224.mae'), ('225.mae'), ('226.mae'), ('227.mae'), ('228.mae'), ('229.mae'), ('230.mae'),
                         #  ('231.mae'), ('232.mae'), ('233.mae'), ('234.mae'), ('235.mae'), ('236.mae'), ('237.mae'), ('238.mae'), ('239.mae'), ('240.mae'),
                         #  ('241.mae'), ('242.mae'), ('243.mae'), ('244.mae'), ('245.mae'), ('246.mae'), ('247.mae'), ('248.mae'), ('249.mae'), ('250.mae'),
                         #  ('251.mae'), ('252.mae'), ('253.mae'), ('254.mae'), ('255.mae'), ('256.mae'), ('257.mae'), ('258.mae'), ('259.mae'), ('260.mae'),
                         #  ('261.mae'), ('262.mae'), ('263.mae'), ('264.mae'), ('265.mae'), ('266.mae'), ('267.mae'), ('268.mae'), ('269.mae'), ('270.mae'),
                         #  ('271.mae'), ('272.mae'), ('273.mae'), ('274.mae'), ('275.mae'), ('276.mae'), ('277.mae'), ('278.mae'), ('279.mae'), ('280.mae'),
                         #  ('281.mae'), ('282.mae'), ('283.mae'), ('284.mae'), ('285.mae'), ('286.mae'), ('287.mae'), ('288.mae'), ('289.mae'), ('290.mae'),
                         #  ('291.mae'), ('292.mae'), ('293.mae'), ('294.mae'), ('295.mae'), ('296.mae'), ('297.mae'), ('298.mae'), ('299.mae'), ('300.mae'),
                          
                         #  ('301.mae'), ('302.mae'), ('303.mae'), ('304.mae'), ('305.mae'), ('306.mae'), ('307.mae'), ('308.mae'), ('309.mae'), ('310.mae'),
                         #  ('311.mae'), ('312.mae'), ('313.mae'), ('314.mae'), ('315.mae'), ('316.mae'), ('317.mae'), ('318.mae'), ('319.mae'), ('320.mae'),
                         #  ('321.mae'), ('322.mae'), ('323.mae'), ('324.mae'), ('325.mae'), ('326.mae'), ('327.mae'), ('328.mae'), ('329.mae'), ('330.mae'),
                         #  ('331.mae'), ('332.mae'), ('333.mae'), ('334.mae'), ('335.mae'), ('336.mae'), ('337.mae'), ('338.mae'), ('339.mae'), ('340.mae'),
                         #  ('341.mae'), ('342.mae'), ('343.mae'), ('344.mae'), ('345.mae'), ('346.mae'), ('347.mae'), ('348.mae'), ('349.mae'), ('350.mae'),
                         #  ('351.mae'), ('352.mae'), ('353.mae'), ('354.mae'), ('355.mae'), ('356.mae'), ('357.mae'), ('358.mae'), ('359.mae'), ('360.mae'),
                         #  ('361.mae'), ('362.mae'), ('363.mae'), ('364.mae'), ('365.mae'), ('366.mae'), ('367.mae'), ('368.mae'), ('369.mae'), ('370.mae'),
                         #  ('371.mae'), ('372.mae'), ('373.mae'), ('374.mae'), ('375.mae'), ('376.mae'), ('377.mae'), ('378.mae'), ('379.mae'), ('380.mae'),
                         #  ('381.mae'), ('382.mae'), ('383.mae'), ('384.mae'), ('385.mae'), ('386.mae'), ('387.mae'), ('388.mae'), ('389.mae'), ('390.mae'),
                         #  ('391.mae'), ('392.mae'), ('393.mae'), ('394.mae'), ('395.mae'), ('396.mae'), ('397.mae'), ('398.mae'), ('399.mae'), ('400.mae'),
                          
                         #  ('401.mae'), ('402.mae'), ('403.mae'), ('404.mae'), ('405.mae'), ('406.mae'), ('407.mae'), ('408.mae'), ('409.mae'), ('410.mae'),
                         #  ('411.mae'), ('412.mae'), ('413.mae'), ('414.mae'), ('415.mae'), ('416.mae'), ('417.mae'), ('418.mae'), ('419.mae'), ('420.mae'),
                         #  ('421.mae'), ('422.mae'), ('423.mae'), ('424.mae'), ('425.mae'), ('426.mae'), ('427.mae'), ('428.mae'), ('429.mae'), ('430.mae'),
                         #  ('431.mae'), ('432.mae'), ('433.mae'), ('434.mae'), ('435.mae'), ('436.mae'), ('437.mae'), ('438.mae'), ('439.mae'), ('440.mae'),
                         #  ('441.mae'), ('442.mae'), ('443.mae'), ('444.mae'), ('445.mae'), ('446.mae'), ('447.mae'), ('448.mae'), ('449.mae'), ('450.mae'),
                         #  ('451.mae'), ('452.mae'), ('453.mae'), ('454.mae'), ('455.mae'), ('456.mae'), ('457.mae'), ('458.mae'), ('459.mae'), ('460.mae'),
                         #  ('461.mae'), ('462.mae'), ('463.mae'), ('464.mae'), ('465.mae'), ('466.mae'), ('467.mae'), ('468.mae'), ('469.mae'), ('470.mae'),
                         #  ('471.mae'), ('472.mae'), ('473.mae'), ('474.mae'), ('475.mae'), ('476.mae'), ('477.mae'), ('478.mae'), ('479.mae'), ('480.mae'),
                         #  ('481.mae'), ('482.mae'), ('483.mae'), ('484.mae'), ('485.mae'), ('486.mae'), ('487.mae'), ('488.mae'), ('489.mae'), ('490.mae'),
                         #  ('491.mae'), ('492.mae'), ('493.mae'), ('494.mae'), ('495.mae'), ('496.mae'), ('497.mae'), ('498.mae'), ('499.mae'), ('500.mae'),                        
                         
                         #  ('501.mae'), ('502.mae'), ('503.mae'), ('504.mae'), ('505.mae'), ('506.mae'), ('507.mae'), ('508.mae'), ('509.mae'), ('510.mae'),
                         #  ('511.mae'), ('512.mae'), ('513.mae'), ('514.mae'), ('515.mae'), ('516.mae'), ('517.mae'), ('518.mae'), ('519.mae'), ('520.mae'),
                         #  ('521.mae'), ('522.mae'), ('523.mae'), ('524.mae'), ('525.mae'), ('526.mae'), ('527.mae'), ('528.mae'), ('529.mae'), ('530.mae'),
                         #  ('531.mae'), ('532.mae'), ('533.mae'), ('534.mae'), ('535.mae'), ('536.mae'), ('537.mae'), ('538.mae'), ('539.mae'), ('540.mae'),
                         #  ('541.mae'), ('542.mae'), ('543.mae'), ('544.mae'), ('545.mae'), ('546.mae'), ('547.mae'), ('548.mae'), ('549.mae'), ('550.mae'),
                         #  ('551.mae'), ('552.mae'), ('553.mae'), ('554.mae'), ('555.mae'), ('556.mae'), ('557.mae'), ('558.mae'), ('559.mae'), ('560.mae'),
                         #  ('561.mae'), ('562.mae'), ('563.mae'), ('564.mae'), ('565.mae'), ('566.mae'), ('567.mae'), ('568.mae'), ('569.mae'), ('570.mae'),
                         #  ('571.mae'), ('572.mae'), ('573.mae'), ('574.mae'), ('575.mae'), ('576.mae'), ('577.mae'), ('578.mae'), ('579.mae'), ('580.mae'),
                         #  ('581.mae'), ('582.mae'), ('583.mae'), ('584.mae'), ('585.mae'), ('586.mae'), ('587.mae'), ('588.mae'), ('589.mae'), ('590.mae'),
                         #  ('591.mae'), ('592.mae'), ('593.mae'), ('594.mae'), ('595.mae'), ('596.mae'), ('597.mae'), ('598.mae'), ('599.mae'), ('600.mae'),

                         #  ('601.mae'), ('602.mae'), ('603.mae'), ('604.mae'), ('605.mae'), ('606.mae'), ('607.mae'), ('608.mae'), ('609.mae'), ('610.mae'),
                         #  ('611.mae'), ('612.mae'), ('613.mae'), ('614.mae'), ('615.mae'), ('616.mae'), ('617.mae'), ('618.mae'), ('619.mae'), ('620.mae'),
                         #  ('621.mae'), ('622.mae'), ('623.mae'), ('624.mae'), ('625.mae'), ('626.mae'), ('627.mae'), ('628.mae'), ('629.mae'), ('630.mae'),
                         #  ('631.mae'), ('632.mae'), ('633.mae'), ('634.mae'), ('635.mae'), ('636.mae'), ('637.mae'), ('638.mae'), ('639.mae'), ('640.mae'),
                         #  ('641.mae'), ('642.mae'), ('643.mae'), ('644.mae'), ('645.mae'), ('646.mae'), ('647.mae'), ('648.mae'), ('649.mae'), ('650.mae'),
                         #  ('651.mae'), ('652.mae'), ('653.mae'), ('654.mae'), ('655.mae'), ('656.mae'), ('657.mae'), ('658.mae'), ('659.mae'), ('660.mae'),
                         #  ('661.mae'), ('662.mae'), ('663.mae'), ('664.mae'), ('665.mae'), ('666.mae'), ('667.mae'), ('668.mae'), ('669.mae'), ('670.mae'),
                         #  ('671.mae'), ('672.mae'), ('673.mae'), ('674.mae'), ('675.mae'), ('676.mae'), ('677.mae'), ('678.mae'), ('679.mae'), ('680.mae'),
                         #  ('681.mae'), ('682.mae'), ('683.mae'), ('684.mae'), ('685.mae'), ('686.mae'), ('687.mae'), ('688.mae'), ('689.mae'), ('690.mae'),
                         #  ('691.mae'), ('692.mae'), ('693.mae'), ('694.mae'), ('695.mae'), ('696.mae'), ('697.mae'), ('698.mae'), ('699.mae'), ('700.mae'),

                         #  ('701.mae'), ('702.mae'), ('703.mae'), ('704.mae'), ('705.mae'), ('706.mae'), ('707.mae'), ('708.mae'), ('709.mae'), ('710.mae'),
                         #  ('711.mae'), ('712.mae'), ('713.mae'), ('714.mae'), ('715.mae'), ('716.mae'), ('717.mae'), ('718.mae'), ('719.mae'), ('720.mae'),
                         #  ('721.mae'), ('722.mae'), ('723.mae'), ('724.mae'), ('725.mae'), ('726.mae'), ('727.mae'), ('728.mae'), ('729.mae'), ('730.mae'),
                         #  ('731.mae'), ('732.mae'), ('733.mae'), ('734.mae'), ('735.mae'), ('736.mae'), ('737.mae'), ('738.mae'), ('739.mae'), ('740.mae'),
                         #  ('741.mae'), ('742.mae'), ('743.mae'), ('744.mae'), ('745.mae'), ('746.mae'), ('747.mae'), ('748.mae'), ('749.mae'), ('750.mae'),                          
                         #  ('751.mae'), ('752.mae'), ('753.mae'), ('754.mae'), ('755.mae'), ('756.mae'), ('757.mae'), ('758.mae'), ('759.mae'), ('760.mae'),
                         #  ('761.mae'), ('762.mae'), ('763.mae'), ('764.mae'), ('765.mae'), ('766.mae'), ('767.mae'), ('768.mae'), ('769.mae'), ('770.mae'),
                         #  ('771.mae'), ('772.mae'), ('773.mae'), ('774.mae'), ('775.mae'), ('776.mae'), ('777.mae'), ('778.mae'), ('779.mae'), ('780.mae'),
                         #  ('781.mae'), ('782.mae'), ('783.mae'), ('784.mae'), ('785.mae'), ('786.mae'), ('787.mae'), ('788.mae'), ('789.mae'), ('790.mae'),
                         #  ('791.mae'), ('792.mae'), ('793.mae'), ('794.mae'), ('795.mae'), ('796.mae'), ('797.mae'), ('798.mae'), ('799.mae'), ('800.mae'),                        
                         
                         #  ('801.mae'), ('802.mae'), ('803.mae'), ('804.mae'), ('805.mae'), ('806.mae'), ('807.mae'), ('808.mae'), ('809.mae'), ('810.mae'),
                         #  ('811.mae'), ('812.mae'), ('813.mae'),
                         # ])
# def test_PlopRotTempGroup(input_file):
#     try:
#         os.remove(os.path.join(REPOS_PATH,'group.dat'))
#         os.remove(os.path.join(REPOS_PATH,"GROUP.dat"))
#     except OSError:
#         pass
#     #res_name = find_resnames_in_mae(os.path.join(TEST_PATH, input_file))[0]
#     subprocess.call([PYTHON_PATH, MAIN_PATH, os.path.join(TEST_PATH, input_file)])
#     subprocess.call([PYTHON_PATH, OLD_MAIN_PATH, os.path.join(TEST_PATH, input_file)])
#     #new_sections = parse_template(os.path.join(REPOS_PATH,res_name.upper()), res_name.upper())
#     #old_sections = parse_template(os.path.join(REPOS_PATH,res_name.lower()), res_name.upper())
    
#     file1 = os.path.join(REPOS_PATH, 'group.dat')
#     file2 = os.path.join(REPOS_PATH, "GROUP.dat")


#     assert filecmp.cmp(file1, file2)

@pytest.mark.parametrize("input_file", [
                          #  (MAE_FILE),
                          # ('MI4.mae'),
                          # ('1o3p_ligand.mae'),
                          ('1.mae'),('2.mae'),('3.mae'),('4.mae'),('5.mae'),('6.mae'),('7.mae'),('8.mae'),('9.mae'),('10.mae'),
                          ('11.mae'),('12.mae'),('13.mae'),('14.mae'),('15.mae'),('16.mae'),('17.mae'),('18.mae'),('19.mae'),('20.mae'),
                          ('21.mae'),('22.mae'),('23.mae'),('24.mae'),('25.mae'),('26.mae'),('27.mae'),('28.mae'),('29.mae'),('30.mae'),
                          # ('31.mae'),('32.mae'),('33.mae'),('34.mae'),('35.mae'),('36.mae'),('37.mae'),('38.mae'),('39.mae'),('40.mae'),
                          # ('41.mae'),('42.mae'),('43.mae'),('44.mae'),('45.mae'),('46.mae'),('47.mae'),('48.mae'),('49.mae'),('50.mae'),
                          # ('51.mae'),('52.mae'),('53.mae'),('54.mae'),('55.mae'),('56.mae'),('57.mae'),('58.mae'),('59.mae'),('60.mae'),
                          # ('61.mae'),('62.mae'),('63.mae'),('64.mae'),('65.mae'),('66.mae'),('67.mae'),('68.mae'),('69.mae'),('70.mae'),
                          # ('71.mae'),('72.mae'),('73.mae'),('74.mae'),('75.mae'),('76.mae'),('77.mae'),('78.mae'),('79.mae'),('80.mae'),
                          # ('81.mae'),('82.mae'),('83.mae'),('84.mae'),('85.mae'),('86.mae'),('87.mae'),('88.mae'),('89.mae'),('90.mae'),
                          # ('91.mae'),('92.mae'),('93.mae'),('94.mae'),('95.mae'),('96.mae'),('97.mae'),('98.mae'),('99.mae'),('100.mae'),

                          # ('100.mae'),('101.mae'),('102.mae'),('103.mae'),('104.mae'),('105.mae'),('106.mae'),('107.mae'),('108.mae'),('109.mae'),('110.mae'),
                          # ('110.mae'),('111.mae'),('112.mae'),('113.mae'),('114.mae'),('115.mae'),('116.mae'),('117.mae'),('118.mae'),('119.mae'),('120.mae'),
                          # ('120.mae'),('121.mae'),('122.mae'),('123.mae'),('124.mae'),('125.mae'),('126.mae'),('127.mae'),('128.mae'),('129.mae'),('130.mae'),
                          # ('131.mae'),('132.mae'),('133.mae'),('134.mae'),('135.mae'),('136.mae'),('137.mae'),('138.mae'),('139.mae'),('140.mae'),
                          # ('141.mae'),('142.mae'),('143.mae'),('144.mae'),('145.mae'),('146.mae'),('147.mae'),('148.mae'),('149.mae'),('150.mae'),
                          # ('151.mae'),('152.mae'),('153.mae'),('154.mae'),('155.mae'),('156.mae'),('157.mae'),('158.mae'),('159.mae'),('160.mae'),
                          # ('161.mae'),('162.mae'),('163.mae'),('164.mae'),('165.mae'),('166.mae'),('167.mae'),('168.mae'),('169.mae'),('170.mae'),
                          # ('171.mae'),('172.mae'),('173.mae'),('174.mae'),('175.mae'),('176.mae'),('177.mae'),('178.mae'),('179.mae'),('180.mae'),
                          # ('181.mae'),('182.mae'),('183.mae'),('184.mae'),('185.mae'),('186.mae'),('187.mae'),('188.mae'),('189.mae'),('190.mae'),
                          # ('191.mae'),('192.mae'),('193.mae'),('194.mae'),('195.mae'),('196.mae'),('197.mae'),('198.mae'),('199.mae'),('200.mae'),

                          # ('201.mae'), ('202.mae'), ('203.mae'), ('204.mae'), ('205.mae'), ('206.mae'), ('207.mae'), ('208.mae'), ('209.mae'), ('210.mae'),
                          # ('211.mae'), ('212.mae'), ('213.mae'), ('214.mae'), ('215.mae'), ('216.mae'), ('217.mae'), ('218.mae'), ('219.mae'), ('220.mae'),
                          # ('221.mae'), ('222.mae'), ('223.mae'), ('224.mae'), ('225.mae'), ('226.mae'), ('227.mae'), ('228.mae'), ('229.mae'), ('230.mae'),
                          # ('231.mae'), ('232.mae'), ('233.mae'), ('234.mae'), ('235.mae'), ('236.mae'), ('237.mae'), ('238.mae'), ('239.mae'), ('240.mae'),
                          # ('241.mae'), ('242.mae'), ('243.mae'), ('244.mae'), ('245.mae'), ('246.mae'), ('247.mae'), ('248.mae'), ('249.mae'), ('250.mae'),
                          # ('251.mae'), ('252.mae'), ('253.mae'), ('254.mae'), ('255.mae'), ('256.mae'), ('257.mae'), ('258.mae'), ('259.mae'), ('260.mae'),
                          # ('261.mae'), ('262.mae'), ('263.mae'), ('264.mae'), ('265.mae'), ('266.mae'), ('267.mae'), ('268.mae'), ('269.mae'), ('270.mae'),
                          # ('271.mae'), ('272.mae'), ('273.mae'), ('274.mae'), ('275.mae'), ('276.mae'), ('277.mae'), ('278.mae'), ('279.mae'), ('280.mae'),
                          # ('281.mae'), ('282.mae'), ('283.mae'), ('284.mae'), ('285.mae'), ('286.mae'), ('287.mae'), ('288.mae'), ('289.mae'), ('290.mae'),
                          # ('291.mae'), ('292.mae'), ('293.mae'), ('294.mae'), ('295.mae'), ('296.mae'), ('297.mae'), ('298.mae'), ('299.mae'), ('300.mae'),
                          
                          # ('301.mae'), ('302.mae'), ('303.mae'), ('304.mae'), ('305.mae'), ('306.mae'), ('307.mae'), ('308.mae'), ('309.mae'), ('310.mae'),
                          # ('311.mae'), ('312.mae'), ('313.mae'), ('314.mae'), ('315.mae'), ('316.mae'), ('317.mae'), ('318.mae'), ('319.mae'), ('320.mae'),
                          # ('321.mae'), ('322.mae'), ('323.mae'), ('324.mae'), ('325.mae'), ('326.mae'), ('327.mae'), ('328.mae'), ('329.mae'), ('330.mae'),
                          # ('331.mae'), ('332.mae'), ('333.mae'), ('334.mae'), ('335.mae'), ('336.mae'), ('337.mae'), ('338.mae'), ('339.mae'), ('340.mae'),
                          # ('341.mae'), ('342.mae'), ('343.mae'), ('344.mae'), ('345.mae'), ('346.mae'), ('347.mae'), ('348.mae'), ('349.mae'), ('350.mae'),
                          # ('351.mae'), ('352.mae'), ('353.mae'), ('354.mae'), ('355.mae'), ('356.mae'), ('357.mae'), ('358.mae'), ('359.mae'), ('360.mae'),
                          # ('361.mae'), ('362.mae'), ('363.mae'), ('364.mae'), ('365.mae'), ('366.mae'), ('367.mae'), ('368.mae'), ('369.mae'), ('370.mae'),
                          # ('371.mae'), ('372.mae'), ('373.mae'), ('374.mae'), ('375.mae'), ('376.mae'), ('377.mae'), ('378.mae'), ('379.mae'), ('380.mae'),
                          # ('381.mae'), ('382.mae'), ('383.mae'), ('384.mae'), ('385.mae'), ('386.mae'), ('387.mae'), ('388.mae'), ('389.mae'), ('390.mae'),
                          # ('391.mae'), ('392.mae'), ('393.mae'), ('394.mae'), ('395.mae'), ('396.mae'), ('397.mae'), ('398.mae'), ('399.mae'), ('400.mae'),
                          
                          # ('401.mae'), ('402.mae'), ('403.mae'), ('404.mae'), ('405.mae'), ('406.mae'), ('407.mae'), ('408.mae'), ('409.mae'), ('410.mae'),
                          # ('411.mae'), ('412.mae'), ('413.mae'), ('414.mae'), ('415.mae'), ('416.mae'), ('417.mae'), ('418.mae'), ('419.mae'), ('420.mae'),
                          # ('421.mae'), ('422.mae'), ('423.mae'), ('424.mae'), ('425.mae'), ('426.mae'), ('427.mae'), ('428.mae'), ('429.mae'), ('430.mae'),
                          # ('431.mae'), ('432.mae'), ('433.mae'), ('434.mae'), ('435.mae'), ('436.mae'), ('437.mae'), ('438.mae'), ('439.mae'), ('440.mae'),
                          # ('441.mae'), ('442.mae'), ('443.mae'), ('444.mae'), ('445.mae'), ('446.mae'), ('447.mae'), ('448.mae'), ('449.mae'), ('450.mae'),
                          # ('451.mae'), ('452.mae'), ('453.mae'), ('454.mae'), ('455.mae'), ('456.mae'), ('457.mae'), ('458.mae'), ('459.mae'), ('460.mae'),
                          # ('461.mae'), ('462.mae'), ('463.mae'), ('464.mae'), ('465.mae'), ('466.mae'), ('467.mae'), ('468.mae'), ('469.mae'), ('470.mae'),
                          # ('471.mae'), ('472.mae'), ('473.mae'), ('474.mae'), ('475.mae'), ('476.mae'), ('477.mae'), ('478.mae'), ('479.mae'), ('480.mae'),
                          # ('481.mae'), ('482.mae'), ('483.mae'), ('484.mae'), ('485.mae'), ('486.mae'), ('487.mae'), ('488.mae'), ('489.mae'), ('490.mae'),
                          # ('491.mae'), ('492.mae'), ('493.mae'), ('494.mae'), ('495.mae'), ('496.mae'), ('497.mae'), ('498.mae'), ('499.mae'), ('500.mae'),                        
                         
                          # ('501.mae'), ('502.mae'), ('503.mae'), ('504.mae'), ('505.mae'), ('506.mae'), ('507.mae'), ('508.mae'), ('509.mae'), ('510.mae'),
                          # ('511.mae'), ('512.mae'), ('513.mae'), ('514.mae'), ('515.mae'), ('516.mae'), ('517.mae'), ('518.mae'), ('519.mae'), ('520.mae'),
                          # ('521.mae'), ('522.mae'), ('523.mae'), ('524.mae'), ('525.mae'), ('526.mae'), ('527.mae'), ('528.mae'), ('529.mae'), ('530.mae'),
                          # ('531.mae'), ('532.mae'), ('533.mae'), ('534.mae'), ('535.mae'), ('536.mae'), ('537.mae'), ('538.mae'), ('539.mae'), ('540.mae'),
                          # ('541.mae'), ('542.mae'), ('543.mae'), ('544.mae'), ('545.mae'), ('546.mae'), ('547.mae'), ('548.mae'), ('549.mae'), ('550.mae'),
                          # ('551.mae'), ('552.mae'), ('553.mae'), ('554.mae'), ('555.mae'), ('556.mae'), ('557.mae'), ('558.mae'), ('559.mae'), ('560.mae'),
                          # ('561.mae'), ('562.mae'), ('563.mae'), ('564.mae'), ('565.mae'), ('566.mae'), ('567.mae'), ('568.mae'), ('569.mae'), ('570.mae'),
                          # ('571.mae'), ('572.mae'), ('573.mae'), ('574.mae'), ('575.mae'), ('576.mae'), ('577.mae'), ('578.mae'), ('579.mae'), ('580.mae'),
                          # ('581.mae'), ('582.mae'), ('583.mae'), ('584.mae'), ('585.mae'), ('586.mae'), ('587.mae'), ('588.mae'), ('589.mae'), ('590.mae'),
                          # ('591.mae'), ('592.mae'), ('593.mae'), ('594.mae'), ('595.mae'), ('596.mae'), ('597.mae'), ('598.mae'), ('599.mae'), ('600.mae'),

                          # ('601.mae'), ('602.mae'), ('603.mae'), ('604.mae'), ('605.mae'), ('606.mae'), ('607.mae'), ('608.mae'), ('609.mae'), ('610.mae'),
                          # ('611.mae'), ('612.mae'), ('613.mae'), ('614.mae'), ('615.mae'), ('616.mae'), ('617.mae'), ('618.mae'), ('619.mae'), ('620.mae'),
                          # ('621.mae'), ('622.mae'), ('623.mae'), ('624.mae'), ('625.mae'), ('626.mae'), ('627.mae'), ('628.mae'), ('629.mae'), ('630.mae'),
                          # ('631.mae'), ('632.mae'), ('633.mae'), ('634.mae'), ('635.mae'), ('636.mae'), ('637.mae'), ('638.mae'), ('639.mae'), ('640.mae'),
                          # ('641.mae'), ('642.mae'), ('643.mae'), ('644.mae'), ('645.mae'), ('646.mae'), ('647.mae'), ('648.mae'), ('649.mae'), ('650.mae'),
                          # ('651.mae'), ('652.mae'), ('653.mae'), ('654.mae'), ('655.mae'), ('656.mae'), ('657.mae'), ('658.mae'), ('659.mae'), ('660.mae'),
                          # ('661.mae'), ('662.mae'), ('663.mae'), ('664.mae'), ('665.mae'), ('666.mae'), ('667.mae'), ('668.mae'), ('669.mae'), ('670.mae'),
                          # ('671.mae'), ('672.mae'), ('673.mae'), ('674.mae'), ('675.mae'), ('676.mae'), ('677.mae'), ('678.mae'), ('679.mae'), ('680.mae'),
                          # ('681.mae'), ('682.mae'), ('683.mae'), ('684.mae'), ('685.mae'), ('686.mae'), ('687.mae'), ('688.mae'), ('689.mae'), ('690.mae'),
                          # ('691.mae'), ('692.mae'), ('693.mae'), ('694.mae'), ('695.mae'), ('696.mae'), ('697.mae'), ('698.mae'), ('699.mae'), ('700.mae'),

                          # ('701.mae'), ('702.mae'), ('703.mae'), ('704.mae'), ('705.mae'), ('706.mae'), ('707.mae'), ('708.mae'), ('709.mae'), ('710.mae'),
                          # ('711.mae'), ('712.mae'), ('713.mae'), ('714.mae'), ('715.mae'), ('716.mae'), ('717.mae'), ('718.mae'), ('719.mae'), ('720.mae'),
                          # ('721.mae'), ('722.mae'), ('723.mae'), ('724.mae'), ('725.mae'), ('726.mae'), ('727.mae'), ('728.mae'), ('729.mae'), ('730.mae'),
                          # ('731.mae'), ('732.mae'), ('733.mae'), ('734.mae'), ('735.mae'), ('736.mae'), ('737.mae'), ('738.mae'), ('739.mae'), ('740.mae'),
                          # ('741.mae'), ('742.mae'), ('743.mae'), ('744.mae'), ('745.mae'), ('746.mae'), ('747.mae'), ('748.mae'), ('749.mae'), ('750.mae'),                          
                          # ('751.mae'), ('752.mae'), ('753.mae'), ('754.mae'), ('755.mae'), ('756.mae'), ('757.mae'), ('758.mae'), ('759.mae'), ('760.mae'),
                          # ('761.mae'), ('762.mae'), ('763.mae'), ('764.mae'), ('765.mae'), ('766.mae'), ('767.mae'), ('768.mae'), ('769.mae'), ('770.mae'),
                          # ('771.mae'), ('772.mae'), ('773.mae'), ('774.mae'), ('775.mae'), ('776.mae'), ('777.mae'), ('778.mae'), ('779.mae'), ('780.mae'),
                          # ('781.mae'), ('782.mae'), ('783.mae'), ('784.mae'), ('785.mae'), ('786.mae'), ('787.mae'), ('788.mae'), ('789.mae'), ('790.mae'),
                          # ('791.mae'), ('792.mae'), ('793.mae'), ('794.mae'), ('795.mae'), ('796.mae'), ('797.mae'), ('798.mae'), ('799.mae'), ('800.mae'),                        
                         
                          # ('801.mae'), ('802.mae'), ('803.mae'), ('804.mae'), ('805.mae'), ('806.mae'), ('807.mae'), ('808.mae'), ('809.mae'), ('810.mae'),
                          # ('811.mae'), ('812.mae'), ('813.mae'),
                         ])
def test_PlopRotTemp_NBNB_Param(input_file):
    try:
      os.remove("lig")
      os.remove("LIG")
    except OSError:
      pass
    # res_name = find_resnames_in_mae(os.path.join(TEST_PATH, input_file))[0]
    res_name = "LIG"
    NAMES = ['sigma', 'epsilon', 'charges', 'SGBr', 'vdwr', 'alpha', 'gamma']
    subprocess.call([PYTHON_PATH, MAIN_PATH, os.path.join(TEST_PATH, input_file)])
    subprocess.call([PYTHON_PATH, OLD_MAIN_PATH, os.path.join(TEST_PATH, input_file)])
    new_sections = parse_template(os.path.join(REPOS_PATH,res_name.upper()), res_name.upper())
    old_sections = parse_template(os.path.join(REPOS_PATH,res_name.lower()), res_name.upper())

    old_names = old_sections[0][3]
    new_names = new_sections[0][3]
  
    NBND_old_sect = old_sections[1]
    NBND_new_sect = new_sections[1]

    NBND_new_ordered = []
    for i in range(len(NBND_old_sect)):
      old_name = old_names[i].strip('_')
      for j, new_params in enumerate(NBND_new_sect):
        if(old_name == new_names[j]):
          NBND_new_ordered.append([param[j] for param in NBND_new_sect])

    for i, (old_section, new_section) in enumerate(zip(zip(*NBND_old_sect), NBND_new_ordered)):
      for old_param, new_param, name in zip(old_section, new_section, NAMES):
        if(old_param != new_param):
          print('[{}, {}, {}, {}]'.format(old_names[i], name, old_param, new_param))
    assert NBND_old_sect == NBND_new_ordered

          



# def test_PlopRotTemp_NBNB_Param(input_file):
#     res_name = input_file.split('.')[0]
#     NAMES = ['sigma', 'epsilon', 'charges', 'SGBr', 'vdwr', 'alpha', 'gamma']
#     subprocess.call([PYTHON_PATH, MAIN_PATH, os.path.join(TEST_PATH, input_file)])
#     subprocess.call([PYTHON_PATH, OLD_MAIN_PATH, os.path.join(TEST_PATH, input_file)])
#     new_sections = parse_template(os.path.join(REPOS_PATH,res_name.upper()), res_name.upper())
#     old_sections = parse_template(os.path.join(REPOS_PATH,res_name.lower()), res_name.upper())



#     old_names = old_sections[0][3]
#     new_names = new_sections[0][3]
  
#     NBND_old_sect = old_sections[1]
#     NBND_new_sect = new_sections[1]

#     for i in zip(*arr)]

#     NBND_new_ordered = []
#     for i, params in enumerate(NBND_old_sect):
#       old_name = old_names[i].strip('_')
#       for j, new_params in enumerate(NBND_new_sect):
#         if(old_name == new_names[j]):
#           NBND_new_ordered.append([params[j] for params in NBND_new_sect])
#           break

#     NBND_old_sect_ordered = []
#     for j in range(len(NBND_old_sect)):
#         NBND_old_sect_ordered.append([params[j] for params in NBND_old_sect])

#     for old_section, new_section in zip(NBND_old_sect_ordered, NBND_new_ordered):
#       for i, (old_param, new_param, name) in enumerate(zip(old_section, new_section, NAMES)):
#         if(old_param != new_param):
#           print('[{}, {}, {}, {}]'.format(old_names[i], name, old_param, new_param))
#     assert NBND_old_sect == NBND_new_ordered
