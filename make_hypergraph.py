from halp.directed_hypergraph import DirectedHypergraph
import gpr_mapping as GPR
from halp.utilities.directed_graph_transformations import to_networkx_digraph
import networkx as nx
from igraph import *
import numpy as np

H = DirectedHypergraph()

met_map = GPR.get_metabolite_associations('./RECON1.json')
for reaction in met_map:
	source = met_map[reaction][0]
	target = met_map[reaction][1]
	H.add_hyperedge(set(source), set(target), weight=1)

g = to_networkx_digraph(H)

#graph = ig.Graph(len(g), zip(*zip(*nx.to_edgelist(g))[:2]))
# print nx.to_edgelist(g)
graph = igraph.Graph.Adjacency((nx.to_numpy_matrix(g) > 0).tolist())
print (graph)
print(graph.vs)
vNames = ["4hpro_LT_m","2pglyc_c","2aobut_m","glyclt_x","4hpro_LT_r","3mox4hpac_c","tststeroneglc_c","Rtotal2coa_c","pyr_m","tststeroneglc_e","pyr_c","pyr_e","Rtotal2coa_m","cs_d_pre4_g","pyr_x","xol7a_r","cmpacna_g","4mptnl_m","cmpacna_c","3dphb_m","dmhptcrn_e","retnglc_r","ala__L_m","ala__L_l","cs_c_pre3_g","dmhptcrn_c","dmhptcrn_m","octa_e","ala__L_e","ksii_core4_deg2_l","ala__L_c","pnto__R_m","octa_c","retnglc_c","L2aadp6sa_m","met__L_c","ala__L_x","cs_d_pre2_g","1mpyr_c","L2aadp6sa_c","ksi_deg2_l","hs_l","uri_n","adn_e","pram_c","ksii_core2_g","galgalfucfucgalacglcgalacglcgal14acglcgalgluside_hs_c","ksi_pre26_g","adn_c","m8mpdol__L_r","cs_e_pre3_g","digalsgalside_hs_c","tdchola_x","digalsgalside_hs_e","3ityr__L_c","atp_e","ump_c","dcholcoa_x","ump_g","atp_c","ump_e","atp_m","atp_n","ump_n","forcoa_x","m7masnB_g","ump_m","ump_r","sgalside_hs_c","forcoa_c","atp_r","m7masnB_r","sgalside_hs_l","atp_x","ksi_pre24_g","omeprazole_e","3hbcoa_m","cbl2_c","pristanal_x","estriolglc_c","dedolp__L_c","cbl2_m","estriolglc_e","m8masn_g","3hbcoa_x","gua_c","gpail_hs_c","din_e","dna_n","glcur1p_c","pristanal_c","R1coa_hs_r","acgagbside_hs_c","acgagbside_hs_g","R2coa_hs_r","acgagbside_hs_l","fucgal14acglcgalgluside_hs_g","urcan_c","gtp_e","gtp_c","cyst__L_c","gtp_n","gtp_m","fuc12gal14acglcgalgluside_hs_g","vitd3_m","mn_l","ptdcacoa_m","3mop_m","vitd3_e","ptdcacoa_c","vitd3_c","dmpp_x","leuktrB4_r","prgstrn_c","cmpacna_n","nmthsrtn_c","prgstrn_r","34hpp_m","cbl1_m","aqcobal_e","selcys_c","aqcobal_c","ncam_c","dhlpro_m","dheas_r","vitd2_c","acglcgalgluside_hs_g","h_x","dheas_e","Nmelys_n","dheas_c","f5hoxkyn_c","cbp_c","asn__L_m","ahandrostan_r","xoltriol_r","asn__L_c","cbp_m","cbp_r","ksi_deg38_l","xoltriol_m","xoltriol_c","glac_c","dhcholestancoa_x","nrvnccoa_c","5aizc_c","nrvnccoa_m","m3mpdol_U_c","h2o2_x","dhcholestancoa_r","i_c","h2o2_c","h2o2_e","ksii_core2_deg6_l","h2o2_n","nrvnccoa_x","h2o2_l","h2o2_m","sl__L_c","acn23acngalgbside_hs_g","acn23acngalgbside_hs_e","leuktrB4_c","acn23acngalgbside_hs_c","sl__L_e","pmtcrn_c","eicostet_c","ksii_core2_deg4_l","cholp_c","cholp_g","clpnd_e","clpnd_c","ksii_core2_deg2_l","cholp_l","cysam_c","acglu_m","mercppyr_m","3h26dm5coa_m","thmtp_e","ksi_deg7_l","thmtp_c","xtsn_c","ttdca_c","met__L_e","glyc_e","glyc_c","ura_e","elaid_c","glyc_m","elaid_e","retinal_cis_13_c","retnglc_e","5hoxindact_c","malcoa_m","tre_c","tre_e","cs_c_deg2_l","core2_g","dopasf_e","m6mpdol_U_r","dopasf_c","chsterol_g","chsterol_e","chsterol_c","m4mpdol__L_c","3php_c","chsterol_m","chsterol_l","ser__L_c","ser__L_e","chsterol_r","m4mpdol__L_r","ac_c","ac_e","m3gacpail_hs_r","ac_g","gd1c_hs_c","lnlncgcrn_c","gd1c_hs_g","gd1c_hs_e","gq1balpha_hs_g","gq1balpha_hs_e","gq1balpha_hs_c","acetone_e","tcynt_m","hs_pre3_g","7dhf_l","7dhf_m","acetone_c","gthox_e","acetone_m","tcynt_e","7dhf_c","7dhf_e","acngalacglcgalgluside_hs_g","2425dhvitd3_m","thrnt_c","acald_m","acald_c","2425dhvitd3_c","crmp_hs_e","acald_x","3uib_c","crmp_hs_c","acald_r","apoC_m","xmp_c","3ohdcoa_x","ksi_deg35_l","prgnlones_r","succ_r","prgnlones_c","cs_b_pre5_g","pre_prot_r","itaccoa_m","nh4_x","nh4_r","core8_g","nh4_n","ins_e","nh4_m","minohp_n","ins_c","ins_l","ins_m","nh4_e","nh4_c","cyan_e","cs_b_deg2_l","cs_c_deg4_l","strch2_e","10fthf6glu_m","2dr5p_c","core6_l","core6_g","mi1456p_n","dgcholcoa_x","tethex3coa_c","prostge2_r","retinal_cis_9_c","cs_a_b_e_pre1_g","tethex3coa_x","prostge2_e","m5masnC_g","prostge2_c","gdp_g","gdp_e","gdp_c","cca_d3_m","gdp_n","gdp_m","cca_d3_c","cca_d3_e","tdchola_e","5aop_c","tetpent3coa_m","tetpent3coa_c","5mta_c","5aop_m","ksi_deg33_l","dad_2_l","ddcacoa_c","dad_2_c","12ppd__S_c","dad_2_e","cspg_a_l","acnacngal14acglcgalgluside_hs_c","tdchola_c","acnacngal14acglcgalgluside_hs_g","tettet6coa_m","peplys_n","antipyrene_c","antipyrene_e","ksi_deg37_l","4fumacac_c","5hxkyn_c","3ohxccoa_x","ksi_deg31_l","cholate_e","digalsgalside_hs_g","13_cis_retnglc_c","retncoa_c","13_cis_retnglc_e","cbasp_c","ksi_pre15_g","emem2gacpail_hs_r","13_cis_retnglc_r","ksi_pre17_g","csn_c","csn_e","etfrd_m","bildglcur_r","3mob_m","ptdca_c","galfucgalacglcgal14acglcgalgluside_hs_e","galfucgalacglcgal14acglcgalgluside_hs_c","3mob_c","htaxol_c","bildglcur_c","gal1p_c","bildglcur_e","5fthf_c","5fthf_e","tetpent3crn_m","prist_x","Lcyst_c","pro__D_x","m8mpdol_U_r","dolmanp_U_r","pro__D_c","pro__D_e","prist_c","dolmanp_U_c","oagt3_hs_g","creat_c","xan_c","dgsn_l","leuktrB4wcooh_r","hnifedipine_c","ksi_deg41_l","creat_m","xan_x","didp_c","hxan_l","galthcrm_hs_g","ksi_pre22_g","hxan_e","didp_n","hxan_c","didp_m","retinal_c","Lcyst_m","hxan_x","spmd_c","prostge1_c","5moxact_c","galacglc13galacglcgal14acglcgalgluside_hs_g","prostge1_e","ppa_c","ppa_e","ppa_m","lac__L_e","bz_e","dcmp_n","dcmp_m","dcmp_l","dcmp_c","sgalside_hs_g","L_dpchrm_c","gbside_hs_c","gbside_hs_e","gbside_hs_g","gbside_hs_l","fuc13galacglcgal14acglcgalgluside_hs_c","malACP_m","ppi_m","ksi_pre11_g","dolichol__L_c","4nphsf_e","ppi_c","malACP_c","fuc__L_l","cs_b_l","4nphsf_c","ksi_pre13_g","ppi_x","acnacngalgbside_hs_e","acnacngalgbside_hs_c","dolichol__L_r","Rtotal3crn_c","ppi_r","pheacgln_e","fuc13galacglcgal14acglcgalgluside_hs_e","pheacgln_c","fprica_c","xu1p__D_c","5hoxnfkyn_c","1p3h5c_m","cholcoas_r","idp_c","cholcoas_x","adrnl_e","m8masn_r","glu5sa_c","andrstrn_c","gal_l","meoh_r","estriol_r","fucfucgalacglc13galacglcgal14acglcgalgluside_hs_g","gal_c","gal_e","glu5sa_m","mag_hs_c","andrstrn_r","urea_c","meoh_c","meoh_l","gar_c","tetpent6coa_m","ppmi1346p_c","ppmi1346p_n","tetpent6coa_c","cytd_c","tetpent6coa_x","cytd_e","cytd_n","cytd_m","cytd_l","udpxyl_c","udpxyl_g","m2emgacpail_hs_r","occoa_x","ksii_core2_deg9_l","pgp_hs_c","udpxyl_r","6a2ohxnt_x","gua_e","12harachd_r","cs_e_l","9_cis_retfa_e","3hibutcoa_m","icit_c","em3gacpail_hs_r","thm_c","estriolglc_r","3snpyr_m","glutcoa_m","galgluside_hs_g","icit_x","trdox_m","fru_e","trdox_c","galgluside_hs_l","crtn_c","hpdcacrn_m","hpdcacrn_c","acg5sa_m","lnlncacrn_m","3dsphgn_c","3padsel_c","thm_m","adpglc_c","succoa_m","m2gacpail_hs_r","adn_m","adn_l","oh1_m","dmantipyrine_c","oh1_c","pnto__R_e","melatn_c","pnto__R_c","dmantipyrine_e","oh1_e","fucfuc12gal14acglcgalgluside_hs_e","fucfuc12gal14acglcgalgluside_hs_g","fucfuc12gal14acglcgalgluside_hs_c","ksi_pre18_g","h_r","cs_pre_g","lcts_l","lcts_c","lcts_g","lcts_e","hretn_c","h_g","h_e","h_c","hpdcacoa_c","h_n","chlstol_r","h_l","h_m","hretn_n","6thf_c","malcoa_c","dolp__L_r","6thf_e","adp_r","thcholstoic_x","adp_x","6thf_m","6thf_l","adp_e","gthox_m","dtdp4d6dm_c","dolp__L_c","adp_c","malcoa_x","thcholstoic_m","adp_n","gthox_c","digalside_hs_g","digalside_hs_c","digalside_hs_l","dmgly_c","fucacngalacglcgalgluside_hs_e","fucacngalacglcgalgluside_hs_g","adprbp_e","adprbp_c","fucacngalacglcgalgluside_hs_c","pheme_m","pheme_e","tdeACP_c","pheme_c","cholate_c","dhocholoylcoa_x","xol7aone_r","lneldccrn_m","fum_c","2mop_m","lxser_g","lneldccrn_c","fum_m","arachd_r","2ameph_c","r5p_m","arachd_e","r5p_c","arachd_c","coucoa_m","carn_c","phllqne_c","phllqne_e","xtp_c","hspg_l","ade_l","hspg_g","hspg_e","ade_c","ade_e","cs_d_l","galt_c","odecoa_c","ppp9_m","leuktrE4_c","whttdca_e","whttdca_c","odecoa_m","nrvnccrn_m","acgalfucgalacglcgalgluside_hs_g","bhb_c","18harachd_r","bhb_e","air_c","odecoa_x","bhb_m","nrvnccrn_c","ocdcyaACP_c","m4masn_g","m3emgacpail_hs_r","ksi_pre3_g","ile__L_c","gthrd_r","cs_d_deg6_l","ocdca_e","gthrd_m","hyptaur_c","ocdca_c","dump_n","gthrd_e","gthrd_c","3aib__D_c","cs_hs_linkage_g","kynate_c","cpppg3_c","selnp_c","crm_hs_l","lac__D_e","thcholst_m","lac__D_c","akgp_hs_c","lac__D_m","estrones_e","crm_hs_c","retinal_11_cis_c","doldp__L_r","tststerones_e","tststerones_c","guln__L_c","phpyr_c","5mthf_e","5mthf_c","ascb__L_c","dhf_c","guln__L_r","retinol_cis_11_e","2maacoa_m","retinol_cis_11_c","acgbgbside_hs_l","retinol_cis_13_c","ptrc_m","ptrc_c","5hoxindoa_c","sphs1p_r","n2m2nmn_l","glu__L_l","34dhoxmand_c","cholcoaone_x","glu__L_c","n2m2nmn_c","glu__L_e","5hoxindoa_m","ficytC_m","glygn3_c","sphs1p_c","sphs1p_e","glu__L_r","betald_m","hmgcoa_r","gt2_hs_g","cholcoads_x","cs_e_deg3_l","ile__L_m","apnnox_e","4mzym_int2_r","hmgcoa_c","cs_e_deg1_l","leuktrB4woh_r","hmgcoa_m","gq1c_hs_g","f26bp_c","ile__L_e","13dpg_c","lneldc_c","hexc_c","hexc_e","lneldc_e","gp1calpha_hs_g","lnlcACP_c","thmmp_e","dgsn_e","thmmp_c","thmmp_m","cs_d_deg2_l","ditp_c","tolbutamide_c","tolbutamide_e","ak2lgchol_hs_c","ditp_n","ak2lgchol_hs_e","ditp_m","ncam_e","msa_m","tsul_m","dudp_m","galfuc12gal14acglcgalgluside_hs_g","trdrd_c","galfuc12gal14acglcgalgluside_hs_e","arg__L_m","trdrd_m","tsul_c","dudp_c","msa_c","cit_m","glcr_c","lnlccoa_c","glcr_m","cit_e","cit_c","g3pc_c","fdp_c","xolest_hs_r","aacoa_x","aacoa_c","xolest_hs_e","hcoumarin_e","xolest_hs_c","hcoumarin_c","aacoa_m","glygn5_e","cs_d_deg5_l","m1mpdol_U_c","m6masnB2_g","25hvitd2_e","25hvitd2_c","25hvitd2_m","aact_m","Rtotal2_c","gd1a_hs_g","dedol_U_c","hdcea_c","l2xser_g","crm_hs_r","ala_B_c","galgalthcrm_hs_g","dolmanp__L_c","tmndnccoa_c","4aabutn_m","dd3coa_m","4aabutn_c","n2m2nmasn_l","urate_c","urate_e","g1m7masnB_r","acgalfucgalacgalfuc12gal14acglcgalgluside_hs_c","urate_x","g1m7masnB_g","aact_c","crm_hs_g","na1_x","xol7ah2_c","whhdca_c","xol7ah2_m","im4act_c","vacccrn_c","xol7ah2_r","vacccrn_m","pecgon_r","taxol_e","nadh_c","asn__L_e","alkylR1oh_x","etoh_e","etoh_c","tetpent3_c","etfox_m","acngalgbside_hs_g","ach_n","4h2oglt_c","glcur_c","5adtststerone_r","acg5p_m","amet_r","tag1p__D_c","3moxtyr_c","hs_deg10_l","5adtststerone_e","5adtststerone_c","glcur_r","amet_m","amet_c","hs_deg16_l","pail34p_hs_c","gt1b_hs_g","mlthf_c","mlthf_m","pail34p_hs_n","strch1_e","5htrp_e","5htrp_c","drib_e","glyc__R_c","4ppan_c","xol7ah_c","dtmp_l","andrstandn_r","ksi_deg13_l","estradiol_c","cs_e_deg4_l","lpam_m","ksii_core4_deg1_l","ksi_deg11_l","drib_c","lyxnt_c","cs_e_deg2_l","estradiol_r","ksi_deg17_l","cmp_e","cmp_g","dcsptn1coa_x","cmp_c","cmp_l","cmp_m","cmp_n","Rtotal3coa_m","ksii_core4_deg3_l","itp_m","cs_e_deg6_l","itp_n","ascb__L_e","i_e","dcsptn1coa_c","17ahprgstrn_c","56dura_c","ksi_deg19_l","17ahprgstrn_r","pppi_n","pppi_m","pppi_c","g1m8mpdol__L_r","pchol_hs_c","man1p_c","pchol_hs_g","pchol_hs_e","g1m6masnB1_g","pchol_hs_m","Rtotal3_c","pchol_hs_r","34dhcinm_c","Rtotal3_e","g1m6masnB1_r","ksii_core2_pre8_g","m5masnB1_g","adhap_hs_c","phaccoa_c","adhap_hs_x","camp_c","camp_e","camp_g","arachdcoa_x","ksii_core2_pre4_g","pdx5p_c","arachdcoa_m","arachdcoa_c","5g2oxpt_x","vacccoa_m","Lpipecol_x","coumarin_e","coumarin_c","vacccoa_c","citmcoa__L_m","hs_deg19_l","ocdcea_c","3mldz_c","1pipdn2c_x","pmtcrn_m","ocdcea_e","apnnox_c","eicostet_e","5thf_l","5thf_m","fmn_c","10fthf5glu_m","10fthf5glu_l","10fthf5glu_c","5thf_e","5thf_c","10fthf5glu_e","fucacngal14acglcgalgluside_hs_g","phom_c","sl__L_m","fe3_e","glu__L_m","acgam_c","hpyr_c","stcrn_m","acgam_e","l2n2m2mn_l","arach_c","arach_e","stcrn_c","acgam_l","m5mpdol__L_r","hpyr_x","fuc132galacglcgal14acglcgalgluside_hs_g","adrn_e","taur_x","adrn_c","cl_e","23doguln_c","cl_c","leuktrD4_r","malttr_c","his__L_c","his__L_e","imp_m","f1p_c","xol7ah2al_m","imp_c","gd3_hs_g","dolp_U_r","dopaqn_c","dolp_U_c","3h26dm5coa_x","gsn_e","13_cis_oretn_n","gsn_c","gsn_l","gsn_m","pylald_c","thym_e","aprut_c","xol24oh_r","ksi_deg40_l","lys__L_n","fucacngal14acglcgalgluside_hs_c","dutp_n","nmn_n","nmn_m","dutp_m","Asn_X_Ser_Thr_l","34dhmald_c","nmn_c","dutp_c","prostg1_c","ttc_ggdp_c","Asn_X_Ser_Thr_r","asp__D_x","tmndnccrn_c","anth_c","m5masnB2_g","pylald_m","tmndnccrn_m","tag__D_c","pyam5p_c","tag__D_e","pyam5p_m","hmgcoa_x","ivcoa_m","hs_deg22_l","hdca_x","xylt_c","tag_hs_c","2obut_c","tag_hs_e","hdca_r","gm3_hs_g","3spyr_m","2mb2coa_m","hdca_e","hdca_c","3spyr_c","hs_deg20_l","ttdca_e","cspg_c_g","10fthf6glu_l","nad_r","3hmp_m","10fthf6glu_c","nad_x","cspg_c_l","10fthf6glu_e","nad_e","cspg_a_e","cspg_a_g","nad_c","nad_m","nm2masn_g","nad_n","ak2g_hs_c","pail5p_hs_n","mmcoa__S_m","ura_c","ksi_pre35_g","mmcoa__S_c","g3m8masn_g","pail5p_hs_c","mmcoa__S_x","cspg_e_e","R4coa_hs_c","g3m8masn_r","Rtotalcrn_c","cspg_e_l","sql_r","pail5p_hs_r","n5m2masn_g","cmpntm2amep_c","ksi_pre31_g","e4p_c","zymstnl_r","gam6p_c","b2coa_m","pcollg5hlys_c","b2coa_x","lnlc_e","lnlc_c","acmana_c","gacpail_hs_c","Tn_antigen_g","val__L_c","gacpail_hs_r","mmcoa__R_m","2dp6mobq_m","abt_c","fucgalacglc13galacglcgal14acglcgalgluside_hs_g","Sfglutth_c","abt_e","fucfucgalacglcgalgluside_hs_e","fucfucgalacglcgalgluside_hs_g","fucfucgalacglcgalgluside_hs_c","Rtotalcrn_m","glc3man_g","pydxn_c","pydxn_e","hs_deg24_l","fucgalacgalfuc12gal14acglcgalgluside_hs_g","chol_r","udpglcur_g","3mgcoa_m","taxol_c","na1_c","na1_e","13_cis_retn_r","na1_g","chol_c","4h2oglt_m","13_cis_retn_n","ach_c","chol_g","ach_e","chol_e","pcreat_m","pcreat_c","micit_c","chol_n","chol_m","13_cis_retn_c","iodine_c","dnad_c","dtmp_c","dd2coa_m","dnad_n","dtmp_m","dnad_m","dadp_n","dadp_m","dadp_c","4tmeabutn_c","tmlys_c","9_cis_retfa_c","adprib_m","1a2425thvitd2_m","tmlys_r","hdd2crn_c","hdd2crn_m","agm_m","1a25dhvitd3_m","xylu__L_c","prpp_c","ipdp_c","ksi_pre34_g","gdchola_c","gdchola_e","val__L_m","gdchola_x","ipdp_x","val__L_e","fucgalfucgalacglcgalgluside_hs_c","fucgalfucgalacglcgalgluside_hs_e","fucgalfucgalacglcgalgluside_hs_g","6hoxmelatn_c","prgnlone_c","accoa_r","prgnlone_m","44mzym_r","accoa_x","g3m8mpdol__L_r","hs_pre15_g","3otdcoa_x","accoa_g","prgnlone_r","accoa_c","accoa_m","glac_r","accoa_n","seasmet_c","hgentis_c","fn2m2masn_g","prostgd2_r","hs_pre13_g","fuc14galacglcgalgluside_hs_c","o2s_x","fuc14galacglcgalgluside_hs_g","fuc14galacglcgalgluside_hs_e","hco3_m","hco3_c","n4m2masn_g","mi34p_c","hco3_e","hs_pre11_g","o2s_n","ser__L_x","o2s_c","alpam_m","o2s_e","ksi_deg20_l","7thf_l","7thf_m","xoltri24_r","7thf_c","7thf_e","xoltri24_e","3htmelys_c","gpi_hs_r","xoltri24_c","hdcea_e","hpdcacoa_m","crtstrn_m","fadh2_r","ksi_deg24_l","fadh2_x","chsterols_c","lpchol_hs_c","fadh2_c","ksi_deg22_l","tdcoa_c","acnam_n","ac_m","acnam_l","chsterols_r","fadh2_m","crtstrn_r","dhcholestanate_m","ksi_deg26_l","estroneglc_r","l2fn2m2masn_g","dhcholestanate_x","zymst_r","paps_c","mpdol_U_c","dhcholestanate_r","ksi_pre32_g","aflatoxin_c","estroneglc_e","ksi_pre30_g","h2co3_m","dudp_n","dtdp_m","t2m26dcoa_x","uppg3_c","retfa_e","dtdp_c","arg__L_e","h2co3_c","retfa_c","gcald_m","arg__L_c","acetol_c","gcald_c","m2mn_c","prostgi2_r","glac_m","m2mn_l","fucgalgbside_hs_e","1mncam_e","cspg_e_g","1mncam_c","pail3p_hs_c","tsul_e","glygn1_c","dna5mtc_n","pail3p_hs_n","dcacoa_c","20ahchsterol_m","gdpddman_c","caro_c","4mtolbutamide_e","o2s_m","4mtolbutamide_c","hista_e","hista_c","pa_hs_m","pa_hs_n","pa_hs_c","caro_e","pa_hs_g","lnlccoa_m","pa_hs_r","pail3p_hs_r","ksi_deg29_l","akdhap_hs_c","gt1alpha_hs_g","akdhap_hs_x","ggn_c","egme_r","limnen_e","cys__L_e","kdn_c","limnen_c","mthgxl_e","mthgxl_c","thr__L_c","tcynt_c","thr__L_e","Nacasp_c","4abut_e","gncore2_g","lnlccrn_m","4abut_c","4abut_l","4abut_m","lnlccrn_c","Nacasp_m","cs_e_deg7_l","2mcit_e","elaidcrn_c","fpram_c","elaidcrn_m","gln__L_e","ksii_core2_pre10_g","gln__L_c","5pmev_x","dag_hs_r","udpg_r","2425dhvitd3_e","orn__D_x","dag_hs_n","udpg_g","dag_hs_c","fucfucfucgalacglc13galacglcgal14acglcgalgluside_hs_c","udpg_c","debrisoquine_e","co2_x","chito2pdol__L_c","adsel_c","co2_r","hs_deg9_l","rbt_c","co2_m","co2_n","co2_c","co2_e","co2_g","eicostetcrn_c","dhdascb_c","dgchol_c","dhdascb_e","dgchol_e","dgchol_x","carveol_c","glx_c","carveol_e","hs_deg7_l","glx_m","s2l2n2m2masn_e","crn_c","ksi_pre6_g","hs_deg1_l","crn_e","dhea_r","hLkynr_c","3amp_l","crn_m","crn_r","ksi_pre8_g","dhea_c","crn_x","hs_deg3_l","mma_c","sel_c","sel_e","tmndnc_e","tmndnc_c","apoC_c","ttccoa_x","aldstrn_m","triodthy_c","mepi_e","aldstrn_c","ksi_pre9_g","amp_c","ksi_pre27_g","cs_e_pre2_g","dgmp_l","crtstrn_c","3mop_c","ocdcaACP_c","dump_m","lald__D_m","whtststerone_r","ksi_pre23_g","tmndnccoa_m","m7masnC_g","11docrtstrn_c","succ_x","ahcys_r","mpdol__L_c","dump_c","Rtotal_l","succ_c","tmndnccoa_x","whtststerone_e","cs_e_pre4_g","ahcys_n","ak2gp_hs_c","5a2opntn_x","11docrtstrn_r","ahcys_c","succ_m","Rtotal_c","acnam_c","man_e","man_g","leuktrC4_r","man_l","Ndmelys_n","leuktrC4_e","ksi_pre1_g","T4hcinnm_m","man_r","cysi__L_e","lgnccoa_x","od2coa_c","nadph_c","nadph_n","nadph_m","od2coa_m","3hxkynam_c","dcdp_n","lgnccoa_m","nadph_r","4hoxpacd_c","k_c","k_e","k_g","dcdp_c","mercplac_c","hs_pre7_g","cmusa_c","xylnt_c","5adtststerones_c","fucgalgbside_hs_g","5adtststerones_e","ksii_core4_pre2_g","prostgf2_c","hs_pre1_g","ksii_core4_pre4_g","galacglcgalacglcgal14acglcgalgluside_hs_g","cs_a_b_pre2_g","oxa_c","fucfuc132galacglcgal14acglcgalgluside_hs_g","fucfuc132galacglcgal14acglcgalgluside_hs_e","fucfuc132galacglcgal14acglcgalgluside_hs_c","34hpl_m","3pg_c","dttOX_c","mercplaccys_e","dttp_n","dttp_m","ps_hs_e","onpthl_c","tyr__L_m","dttp_c","onpthl_e","dhap_x","paps_g","dcyt_n","argsuc_c","but_m","5hoxindact_m","hs_pre8_g","but_e","ksii_core4_pre9_g","retinol_9_cis_e","ap4a_c","but_c","gullac_r","ksii_core2_deg7_l","estroneglc_c","lgnccrn_c","ksii_core2_e","uri_l","uri_m","uri_c","lgnccrn_m","ksii_core2_l","gullac_c","uri_e","hom__L_c","aflatoxin_e","ksii_core2_deg5_l","ksii_core4_l","gm1b_hs_g","ksii_core4_g","ksii_core4_e","ksii_core2_deg3_l","avite1_e","34hpp_c","bdg2hc_c","dtdp_n","c2m26dcoa_x","Rtotal2_e","uacgam_c","adrncrn_c","Rtotalcoa_c","lpro_m","uacgam_g","Rtotalcoa_m","rnam_c","lnlncgcrn_m","adrncrn_m","lald__L_m","cyan_m","kdnp_c","Rtotalcoa_x","dedolp_U_c","thcholoylcoa_x","cyan_c","cs_c_deg1_l","dpcoa_c","n3m2masn_g","mi134p_c","dpcoa_l","palmACP_c","1a25dhvitd2_m","ksi_deg8_l","cs_c_deg3_l","creat_e","gm2_hs_g","nicrnt_n","nicrnt_m","nicrnt_c","emgacpail_hs_r","ibcoa_m","hpdca_e","hpdca_c","decdp_m","estradiolglc_r","gly_l","decdp_c","4izp_c","tettet6crn_c","n2m2mn_l","tettet6crn_m","hs_pre4_g","mi13456p_c","mi13456p_n","mercppyr_c","ala__D_x","ala__D_e","bamppald_m","ala__D_c","tymsf_e","dmgly_m","ksii_core4_pre7_g","leuktrF4_c","ala__D_l","bamppald_c","pap_g","ksi_deg34_l","ksi_pre20_g","pap_c","ksi_deg1_e","ksi_deg1_l","udpgal_c","3aib__D_m","udpgal_g","minohp_c","gm1_hs_g","ksi_deg3_l","3aib__D_e","ksi_deg36_l","ksi_deg5_l","cs_a_deg5_l","estradiolglc_e","adprib_e","adprib_c","estradiolglc_c","acac_m","acac_c","naglc2p__L_c","acac_e","Ser_Gly_Ala_X_Gly_r","cs_a_deg1_l","ksi_pre28_g","pep_c","L2aadp_c","homoval_c","ecgon_r","L2aadp_m","pep_m","5mdru1p_c","im4act_m","cs_a_deg3_l","cs_b_deg1_l","cs_c_deg5_l","omeprazole_c","gq1b_hs_c","11_cis_retfa_c","gq1b_hs_g","11_cis_retfa_e","gq1b_hs_e","lum3_c","em2emgacpail_prot_hs_r","core4_g","hdcoa_c","em2emgacpail_hs_r","7dhchsterol_r","hdcoa_m","glyald_c","urea_m","malt_l","7dhchsterol_c","malt_c","malt_e","glx_x","cs_d_pre5_g","xoldiolone_r","taur_e","nwharg_c","taur_c","s7p_c","xoldiolone_m","xoldiolone_c","peamn_c","lys__L_e","lys__L_c","galside_hs_l","galside_hs_c","lys__L_m","pac_c","galside_hs_g","cs_d_pre3_g","perillyl_c","inost_e","inost_c","cs_c_pre2_g","focytC_m","ksi_pre2_g","dcacoa_x","4abutn_c","fol_e","fol_c","ksi_deg30_l","udpglcur_c","2amac_c","12ppd__R_c","ksi_deg32_l","hs_deg5_l","udpglcur_r","hs_deg18_l","gdpfuc_c","gt1c_hs_g","gdpfuc_g","ser__L_m","duri_m","galfuc12gal14acglcgalgluside_hs_c","duri_n","eaflatoxin_e","for_r","duri_e","prostgh2_c","duri_c","for_e","cs_e_pre5a_g","for_c","for_m","for_n","prostgh2_r","udpacgal_r","prpncoa_x","12RHPET_c","prpncoa_m","udpacgal_g","udpacgal_c","udpacgal_l","ksii_core2_deg8_l","nadp_c","nadp_e","mgacpail_hs_r","nadp_m","nadp_n","galfucgalacglcgalgluside_hs_g","nadp_r","m2mpdol_U_c","nadp_x","mi1346p_n","vitd2_e","dgmp_m","nicrns_c","ksi_g","ksi_e","vitd2_m","cysi__L_c","T4hcinnm_c","dgmp_c","ksi_l","thmpp_m","thmpp_c","dhor__S_c","2mp2coa_m","nadh_x","coa_r","f6p_c","dsmsterol_r","nadh_r","c2m26dcoa_m","coa_x","coa_g","avite1_c","coa_c","nadh_m","hom__L_e","coa_n","coa_m","coa_l","amuco_c","sT_antigen_g","am6sa_c","acrn_m","fucgalacglcgal14acglcgalgluside_hs_g","pail345p_hs_n","fuc1p__L_c","avite2_e","pail345p_hs_c","m2mpdol__L_c","dcsptn1coa_m","mi1p__D_n","mi1p__D_c","ksi_deg39_l","thbpt_c","thbpt_n","oaa_x","m5mpdol_U_r","Lkynr_c","gudac_c","etoh_x","dopa_e","dhcholoylcoa_r","pro__L_l","pro__L_m","dopa_c","dhcholoylcoa_x","pro__L_c","pro__L_e","orot_c","sphmyln_hs_c","oaa_c","alpro_m","oaa_m","pro__L_r","memgacpail_hs_r","fucgalacglcgalgluside_hs_g","ts3_c","arab__L_e","arab__L_c","mn_c","glyb_m","c226coa_c","acACP_c","glyb_e","c226coa_m","glyb_c","m3mpdol__L_c","56dihindlcrbxlt_c","citr__L_m","acglcgalacglcgal14acglcgalgluside_hs_g","thcrm_hs_c","pail4p_hs_r","citr__L_c","thcrm_hs_g","c226coa_x","4tmeabut_m","4tmeabut_c","eandrstrn_r","ocdcya_c","5forthf_c","leuktrD4_c","gtocophe_c","gtocophe_e","nrpphrsf_e","nrpphrsf_c","nrpphr_c","selhcys_c","nrpphr_e","gp1c_hs_g","galacglcgal14acglcgalgluside_hs_g","gp1c_hs_e","adrncoa_x","gp1c_hs_c","eicostetcoa_m","2425dhvitd2_e","sphmyln_hs_l","methf_c","fru_c","2425dhvitd2_c","2425dhvitd2_m","methf_m","eicostetcoa_c","sphmyln_hs_g","Ser_Thr_l","ddcacoa_x","galgalgalthcrm_hs_g","acgal_c","g3m8mpdol_U_r","Ser_Thr_g","selmeth_c","ksi_pre19_g","peplys_e","44mctr_r","forglu_c","pmtcoa_x","35diotyr_c","44mctr_c","pydx5p_c","4hphac_c","4hphac_e","pmtcoa_m","phyt_x","mi4p__D_c","mem2emgacpail_prot_hs_r","pydx_e","2dp6mobq_me_m","pan4p_c","cgly_e","ser__D_c","phyt_c","ser__D_e","phyt_e","cgly_c","core7_g","5adtststeroneglc_r","cs_a_b_pre3_g","3hbcoa__R_m","4aphdob_c","5adtststeroneglc_e","2dr1p_c","5adtststeroneglc_c","andrstrnglc_r","core3_g","5dhf_l","5dhf_m","1glyc_hs_e","sbt__D_c","pi_r","5dhf_e","5dhf_c","cs_b_pre4_g","andrstrnglc_c","core5_g","pi_l","pi_m","pi_n","andrstrnglc_e","pi_c","pi_e","pi_g","whddca_c","whddca_e","grdp_c","hcys__L_c","34dhoxpeg_c","Ntmelys_r","34dhoxpeg_e","hestratriol_c","leuktrA4_c","ctp_n","ctp_m","ctp_c","grdp_x","ksi_pre33_g","pppg9_c","lnlncgcoa_m","lnlncgcoa_c","pppg9_m","xoltetrol_m","lthstrl_r","retinol_9_cis_c","lnlncgcoa_x","itacon_m","ru5p__D_c","glcur_l","retn_n","retn_c","sarcs_e","nm4masn_g","retn_e","3odcoa_x","retn_r","chtn_e","fucgalacgalfucgalacglcgal14acglcgalgluside_hs_g","34dhphe_c","chtn_c","5hxkynam_c","sphgn_c","acgalfucgalacgalfucgalacglcgal14acglcgalgluside_hs_g","acgalfucgalacgalfucgalacglcgal14acglcgalgluside_hs_e","acgalfucgalacgalfucgalacglcgal14acglcgalgluside_hs_c","amet_n","andrstndn_r","sphgn_r","ksi_pre14_g","cspg_d_l","tethex3_e","acglcgalgbside_hs_g","tethex3_c","2dp6mep_m","34dhpac_c","cs_a_l","biocyt_m","cs_e_pre5b_g","cs_d_deg1_l","13dampp_c","biocyt_e","biocyt_c","xylt_e","q10_c","g1m8masn_g","acglcgal14acglcgalgluside_hs_g","gluala_e","g1m8masn_r","fucfucfucgalacglc13galacglcgal14acglcgalgluside_hs_e","fucfucfucgalacglc13galacglcgal14acglcgalgluside_hs_g","oagd3_hs_g","ggdp_c","lipoate_e","oagd3_hs_c","lipoate_c","damp_l","r1p_c","R3coa_hs_c","damp_c","idp_m","idp_n","idp_e","glc__D_r","R2coa_hs_c","naglc2p_U_c","glc__D_l","glc__D_c","3sala_c","glc__D_g","glc__D_e","pydx_c","ebastine_c","frdp_c","ebastine_e","bvite_e","m6masnB1_g","ksii_core2_pre6_g","pecgoncoa_r","ebastine_r","frdp_r","fucfucgalacglcgalacglcgal14acglcgalgluside_hs_g","mi3456p_c","lneldcACP_c","frdp_x","pail34p_hs_r","seln_c","e4hglu_c","gluside_hs_r","thcholstoic_r","akg_x","e4hglu_m","eaflatoxin_c","gluside_hs_g","gluside_hs_c","gluside_hs_l","galacgalfuc12gal14acglcgalgluside_hs_g","pcrn_x","ksi_pre12_g","pcrn_c","Rtotal_x","4ppan_m","mhista_c","pcrn_m","cs_d_deg3_l","cs_c_l","ksi_pre16_g","galfucgalacglcgal14acglcgalgluside_hs_g","ksi_pre10_g","lac__L_m","dolglcp__L_c","pail_hs_r","stcoa_x","lac__L_c","cholcoar_x","stcoa_c","dolglcp__L_r","pail_hs_n","cdpdag_hs_c","cdpdag_hs_m","pail_hs_c","cholcoar_r","ptdca_e","adrnl_c","stcoa_m","estrones_c","gam_c","n2m2nm_l","gam_e","icit_m","etha_c","tetpent3crn_c","pail35p_hs_r","g6p_r","n2m2nm_c","estrones_r","occoa_m","ksii_core2_deg1_l","pail35p_hs_c","g6p_c","occoa_c","id3acald_m","idour_l","3m4hpga_c","htaxol_e","o2_n","phytcoa_c","idour_c","lnlncacoa_m","phytcoa_x","nac_c","lnlncacoa_c","o2_m","4mlacac_c","mal__L_m","ACP_m","g1p_c","mal__L_c","ACP_c","cs_d_deg4_l","11docrtsl_c","oxa_x","11docrtsl_r","dtdp4d6dg_c","gt1a_hs_g","xol27oh_r","ahandrostanglc_c","ahandrostanglc_e","btcoa_m","R6coa_hs_c","xol27oh_m","ahandrostanglc_r","tetpent3_e","Rtotal3crn_m","gmp_c","gmp_e","ethamp_c","gmp_l","gmp_m","gmp_n","debrisoquine_c","ethamp_r","pacald_c","galacglcgalgbside_hs_e","galacglcgalgbside_hs_g","btn_c","hexccrn_m","btn_e","galacglcgalgbside_hs_c","dag_hs_e","hexccrn_c","btn_m","btn_n","cspg_c_e","xol7ah3_m","4hbzcoa_m","56dthm_c","dhcrm_hs_c","thcys_c","4nph_e","mi13p_c","4nph_c","g2m8mpdol__L_r","6dhf_c","acnamp_c","mescoa_m","pro__D_l","nmptrc_c","crvnc_e","crvnc_c","phe__L_e","ca2_e","ca2_c","phe__L_c","phe__L_m","lnlnca_e","tyr__L_e","5mdr1p_c","lnlnca_c","tyr__L_c","3hpp_c","ahdt_n","prostgd2_c","lnlncg_e","lnlncg_c","tym_c","ahdt_c","itp_c","seasmet_n","dmnoncrn_m","rbt_e","cs_c_d_e_pre1_g","dmnoncoa_x","pser__L_c","ppmi12346p_n","ppmi12346p_c","ddcaACP_c","3htmelys_m","gd1b_hs_g","cs_e_deg5_l","dmnoncoa_m","trp__L_e","dmnoncoa_c","acgam1p_c","trp__L_c","sTn_antigen_g","s2l2fn2m2masn_g","s2l2fn2m2masn_e","acgalfucgalacglcgal14acglcgalgluside_hs_g","24nph_c","2mbcoa_m","s2l2fn2m2masn_l","xyl__D_l","perillyl_e","acgal1p_c","pe_hs_m","xyl__D_e","xyl__D_c","dcamp_c","dhf_e","amp_r","tststeroneglc_r","dhf_l","dhf_m","amp_x","amp_g","orn_e","amp_e","orn_c","ha_pre1_e","mi145p_c","orn_m","amp_m","ha_pre1_l","mescon_m","aps_c","xoltri25_e","35cgmp_n","dcyt_c","eicostetcrn_m","tetpent6_c","tetpent6_e","dcyt_e","acald_e","dcyt_m","dcyt_l","sphings_r","fucfucfucgalacglcgal14acglcgalgluside_hs_c","fucfucfucgalacglcgal14acglcgalgluside_hs_g","ntm2amep_c","fucfucfucgalacglcgal14acglcgalgluside_hs_e","clpn_hs_c","sphings_c","apoC_Lys_c","apoC_Lys_m","2hyoxplac_c","sphings_l","3oddcoa_x","so3_c","ind3ac_c","so3_m","s2l2n2m2masn_l","dha_c","biliverd_c","dedol__L_c","10fthf_e","fald_l","fald_m","dlnlcgcoa_c","dlnlcgcoa_m","xser_r","fald_c","fald_x","mi1345p_c","mi1345p_n","xser_g","pmtcoa_c","vacc_e","fad_x","mi1346p_c","mercplaccys_c","hs_deg14_l","fad_r","fad_m","man6p_c","hs_deg23_l","fad_c","apoC_Lys_btn_m","ebastineoh_r","ebastineoh_e","ksi_deg15_l","ebastineoh_c","glygn2_c","glygn2_e","xoldioloneh_c","so4_r","gd3_hs_c","acgal_l","hestratriol_r","tettet6coa_c","acgal_g","cs_a_deg4_l","so4_c","tettet6coa_x","hmbil_c","4mzym_int1_r","so4_e","glucys_c","urea_e","so4_l","so4_m","mi3p__D_c","thyox__L_c","prostgd2_e","thyox__L_e","octeACP_c","23dpg_c","Ssq23epx_r","4mop_m","c226crn_m","4mop_c","c226crn_c","Ntmelys_n","acngalacglcgal14acglcgalgluside_hs_e","acngalacglcgal14acglcgalgluside_hs_g","sarcs_m","melanin_c","acngalacglcgal14acglcgalgluside_hs_c","bilirub_e","hexccoa_m","dtt_c","sarcs_c","bilirub_c","glu5p_m","gm1a_hs_g","sarcs_x","bilirub_r","o2_c","utp_c","utp_e","o2_e","dcsptn1crn_c","utp_m","hdcecrn_m","utp_n","o2_r","g2m8mpdol_U_r","arachdcrn_m","gt1a_hs_e","48dhoxquin_c","o2_x","hxdcal_r","gt1a_hs_c","dem2emgacpail_prot_hs_r","dak2gpe_hs_c","dmnoncrn_x","ala_B_e","6dhf_e","xol7aone_c","ala_B_m","6dhf_m","6dhf_l","txa2_e","txa2_c","hdcecrn_c","arachdcrn_c","ind3ac_m","txa2_r","f1a_g","35cgmp_e","bz_r","35cgmp_g","f1a_l","35cgmp_c","saccrp__L_m","hexccoa_x","oagd3_hs_e","glygn4_e","clpndcrn_c","ps_hs_g","hexccoa_c","25hvitd3_m","arg__D_x","zym_int2_r","retinol_c","25hvitd3_e","retinol_e","25hvitd3_c","clpndcrn_m","h2o_c","3hpcoa_x","h2o_e","tststerone_c","h2o_g","ps_hs_c","nac_e","h2o_m","h2o_l","h2o_n","h2o_r","tststerone_r","h2o_x","spc_hs_c","spc_hs_e","bz_c","triodthy_e","acac_x","apoC_Lys_btn_c","cspg_d_g","malttr_l","m_em_3gacpail_hs_r","pail45p_hs_n","cspg_b_l","malttr_e","pail45p_hs_c","l2n2m2mn_c","cspg_b_g","cspg_b_e","sph1p_r","g1m7masnC_r","2pg_c","pcollglys_c","sph1p_c","sph1p_e","triodthysuf_e","ttdcea_c","3mox4hoxm_c","glyc3p_c","glyc3p_m","aldstrn_e","10fthf_l","10fthf_m","1pyr5c_m","1pyr5c_c","normete__L_c","mepi_c","10fthf_c","avite2_c","cspg_d_e","17ahprgnlone_c","xu5p__D_c","gd1alpha_hs_g","3dpdhb_m","ribflv_e","17ahprgnlone_r","12HPET_c","ribflv_c","11docrtsl_m","m4mpdol_U_r","leuktrA4_r","gln__L_m","m4mpdol_U_c","dtdprmn_c","pi_x","25aics_c","ksii_core2_pre2_g","fuc__L_c","m_em_3gacpail_prot_hs_r","wharachd_r","fuc__L_e","3mlda_e","3mlda_c","eryth_c","ha_deg1_l","dkmpp_c","phpyr_m","ppi_n","5oxpro_c","cdpea_c","4pyrdx_e","glyc__S_e","glyc__S_c","ksii_core4_pre10_g","hs_deg13_l","gdpmann_c","lgt__S_m","n2m2masn_g","ksi_deg12_l","lgt__S_c","pe_hs_r","m6mpdol__L_r","ksi_deg10_l","myrsACP_c","dhbpt_c","pe_hs_e","pe_hs_g","pe_hs_c","acnacngalgbside_hs_g","ps_hs_r","34dhpha_c","3ohodcoa_x","btamp_c","thf_e","thf_c","ksi_deg14_l","ga2_hs_g","thf_l","thf_m","btamp_m","clpndcoa_m","adp_m","m6masnC_g","ksii_core4_deg4_l","1a2425thvitd3_m","ksi_deg18_l","clpndcoa_c","fucgal14acglcgalgluside_hs_e","m6masnA_g","xolest2_hs_c","lneldccoa_c","fucgal14acglcgalgluside_hs_c","clpndcoa_x","strdnccrn_c","ksi_deg16_l","lneldccoa_m","xolest2_hs_e","gm2a_hs_g","mev__R_r","11docrtstrn_m","g1m7masnC_g","bilglcur_e","mev__R_x","pristcoa_x","tettet6_c","tettet6_e","bilglcur_c","pristcoa_c","hdeACP_c","bilglcur_r","nrvnc_c","nrvnc_e","fuc13galacglcgal14acglcgalgluside_hs_g","tchola_e","glc1man_g","thym_c","acnacngal14acglcgalgluside_hs_e","ksii_core2_pre9_g","acmanap_c","gpi_sig_r","ak2gpe_hs_c","galacglcgalgluside_hs_g","pglyc_hs_c","pglyc_hs_e","adrncoa_m","adrncoa_c","ksii_core2_pre7_g","dedoldp_U_c","whtststerone_c","24nph_e","no_e","dimp_c","T_antigen_g","no_c","succ_e","rib__D_e","ksii_core2_pre5_g","rib__D_c","ahcys_m","dgdp_c","dedoldp__L_c","ksi_pre25_g","dgdp_m","dgdp_n","ksii_core2_pre1_g","Rtotal_e","oagt3_hs_e","oagt3_hs_c","hs_deg17_l","peracd_c","Lfmkynr_c","chito2pdol_U_c","peracd_m","m7mpdol__L_r","id3acald_c","appnn_c","tetpent6crn_c","thm_e","fe2_m","fe2_c","tetpent6crn_m","fe2_e","acgalfuc12gal14acglcgalgluside_hs_g","appnn_e","man_c","hs_deg15_l","m7mpdol_U_r","co_c","leu__L_c","leu__L_e","n4abutn_c","co_e","adocbl_m","n4abutn_m","leu__L_m","thbpt4acam_c","3aib_c","betald_c","3aib_e","andrstrn_e","acgbgbside_hs_g","acgbgbside_hs_c","3aib_m","dgsn_c","gp1calpha_hs_c","ppbng_c","5homeprazole_c","cdp_n","cdp_m","gp1calpha_hs_e","cdp_c","dgsn_m","fucacngal14acglcgalgluside_hs_e","n2m2nmasn_g","n2m2nmasn_e","dolichol_U_c","leuktrC4_c","pail4p_hs_c","acgalfucgalacgalfuc12gal14acglcgalgluside_hs_g","pail4p_hs_n","acgalfucgalacgalfuc12gal14acglcgalgluside_hs_e","Tn_antigen_l","g3p_c","paf_hs_c","dolichol_U_r","dolmanp__L_r","paf_hs_e","ptth_c","tetpent3coa_x","mem2emgacpail_hs_r","g2m8masn_r","ga1_hs_g","acgpail_hs_c","meoh_e","gchola_c","gchola_e","gchola_x","4hdebrisoquine_c","4hdebrisoquine_e","g2m8masn_g","Tyr_ggn_c","selcyst_c","Rtotal3coa_c","3hmbcoa_m","thymd_l","thymd_m","xoltri27_c","thymd_e","mag_hs_e","thymd_c","34dhphe_e","estriol_c","46dhoxquin_c","udp_r","Ser_Gly_Ala_X_Gly_l","npthl_c","npthl_e","ump_l","udp_c","udp_e","udp_g","udp_m","udp_l","udp_n","ksi_deg28_l","fgam_c","aicar_c","lgnc_c","lgnc_e","hs_deg21_l","dcsptn1_e","dcsptn1_c","dcdp_m","estrone_c","quln_c","s2l2n2m2mn_l","estrone_r","5dpmev_x","s2l2n2m2mn_c","adpman_c","glc2man_g","tdcoa_x","orot5p_c","tdcoa_m","acn13acngalgbside_hs_c","lpchol_hs_e","acn13acngalgbside_hs_e","lgnccoa_c","acn13acngalgbside_hs_g","2hb_c","2hb_e","R1coa_hs_c","3mb2coa_m","2coum_c","strdnccrn_m","t2m26dcoa_m","fna5moxam_c","ksi_pre36_g","xylu__D_c","nadph_x","hs_deg6_l","cys__L_c","3dpdhb_me_m","ksii_core2_pre3_g","15HPET_c","2mcit_c","ppcoa_x","gmp_g","cys__L_m","hs_deg8_l","aprgstrn_c","aprgstrn_e","ppcoa_c","oretn_c","ppcoa_m","oretn_n","hpyr_m","ak2gchol_hs_c","6pgl_c","hs_deg2_l","gly_x","acngal14acglcgalgluside_hs_g","thp2c_x","thp2c_c","6pgl_r","gly_m","phyt2ohcoa_x","hs_deg4_l","prostgf2_e","gly_e","gly_c","lald__D_c","nformanth_c","acgam6p_c","2mcacn_c","dhap_c","ksi_pre5_g","oxa_e","rbl__D_c","sprm_c","nifedipine_e","nifedipine_c","ksi_pre7_g","fucgalgbside_hs_c","dsT_antigen_g","sucsal_m","lald__L_c","6pgc_c","ksii_core4_pre3_g","lanost_c","lnlncacrn_c","hs_deg25_l","4hbz_m","tymsf_c","6pgc_r","lanost_r","arachcoa_x","gal14acglcgalgluside_hs_g","dgpi_prot_hs_r","arachcoa_m","arachcoa_c","asp__D_c","13_cis_oretn_c","asp__D_e","whhdca_e","din_c","dgtp_c","ksi_pre4_g","cortsn_r","5homeprazole_e","dgtp_m","dgtp_n","acrn_c","dolglcp_U_r","gd1b2_hs_e","gd1b2_hs_g","gd1b2_hs_c","dolglcp_U_c","acrn_r","acrn_x","galgalgalthcrm_hs_c","gpi_prot_hs_r","asp__L_m","asp__L_c","asp__L_e","pepslys_r","Rtotal2crn_m","ksi_deg4_l","adrn_x","cs_b_deg3_l","Rtotal2crn_c","cdpchol_c","srtn_c","strdnc_e","strdnc_c","srtn_e","pydx5p_m","hs_deg11_l","4ppcys_c","hs_pre6_g","alpa_hs_c","ksii_core4_pre5_g","alpa_hs_m","dmhptcoa_m","bvite_c","lys__L_x","datp_n","hs_pre14_g","datp_m","datp_c","R5coa_hs_c","dmhptcoa_c","3sala_m","sucr_e","q10_m","dcsptn1crn_m","im4ac_c","im4ac_m","hs_pre10_g","galgalfucfucgalacglcgalacglcgal14acglcgalgluside_hs_g","galgalfucfucgalacglcgalacglcgal14acglcgalgluside_hs_e","amp_l","seahcys_c","hs_pre12_g","3dhguln_c","seahcys_n","dmpp_c","biocyt_n","3hanthrn_c","ptdcacrn_m","2kmb_c","Nacsertn_c","5HPET_c","acglc13galacglcgal14acglcgalgluside_hs_g","ptdcacrn_c","dtdpglu_c","galgbside_hs_g","3snpyr_c","xoltri25_r","galgalgalthcrm_hs_e","2dpmhobq_m","hs_pre9_g","xoltri25_c","dlnlcg_c","g1m8mpdol_U_r","gncore1_g","Tyr_ggn_e","gd2_hs_g","mi145p_n","coke_r","ksi_deg23_l","ksi_deg21_l","vacc_c","doldp_U_r","ksi_deg27_l","4pyrdx_c","q10h2_c","pydam_c","galgluside_hs_c","q10h2_m","dlnlcg_e","ddsmsterol_r","pydam_e","ksi_deg25_l","42A3HP24DB_c","ps_hs_m","imp_e","ddca_c","hs_pre2_g","2c23dh56dhoxin_c","1p2cbxl_x","3hpcoa_c","ksii_core4_pre1_g","6pthp_c","6pthp_n","triodthysuf_c","ru5p__D_r","3hpcoa_m","mi14p_c","mi14p_n","ksi_deg9_l","ametam_c","glyclt_c","6htststerone_c","tchola_c","glyclt_m","6htststerone_e","tchola_x","hestratriol_e","6htststerone_r","akg_c","cmp2amep_c","dlnlcgcrn_m","akg_e","cs_a_deg2_l","dlnlcgcrn_c","akg_m","pd3_c","akg_r","ha_l","cholcoa_x","ha_e","ksii_core4_pre6_g","xoltri27_e","gt3_hs_c","4hglusa_m","hs_deg12_l","gt3_hs_g","xoltri27_r","ksii_core4_pre8_g","xol25oh_r","s2l2n2m2m_l","s2l2n2m2m_c","hs_pre5_g","10fthf7glu_c","10fthf7glu_e","ksi_deg6_l","10fthf7glu_m","10fthf7glu_l","arachcrn_c","fucfucgalacglcgal14acglcgalgluside_hs_g","arachcrn_m","m7masnA_g","ksi_pre21_g","cala_c","fucacgalfucgalacglcgalgluside_hs_g","fucacgalfucgalacglcgalgluside_hs_e","fucacgalfucgalacglcgalgluside_hs_c","m1mpdol__L_c","strdnccoa_c","crtsl_c","ksi_pre29_g","strdnccoa_m","crtsl_m","crtsl_r","galacgalfucgalacglcgal14acglcgalgluside_hs_g","dhlam_m","acorn_c","trypta_c","strdnccoa_x","dctp_m","dctp_n","2oxoadp_c","hdd2coa_x","2oxoadp_m","dctp_c","hdd2coa_m","ttdcrn_c","odecrn_m","dxtrn_c","odecrn_c","hdd2coa_c","ttdcrn_m"]
adjMatrix = (nx.to_numpy_matrix(g) > 0).tolist()
graph = ig.Graph.Adjacency(adjMatrix)
graph.vs['name'] = vNames
print graph

#USED TO GENERATE vNames 
#name_dict = g.nodes(data=True)
# vertex_names = set()
# for e in graph.es: 
# 	s = name_dict[e.source][0]
# 	t = name_dict[e.target][0]
# 	vertex_names.update([s,t])
# vertices = list(vertex_names)
# f = open('dump.txt','w')
# for x in vertices: 
# 	f.write("\"" + str(x)+ "\",")
# f.close()
# # print vertices 
# graph.colnames = vertices