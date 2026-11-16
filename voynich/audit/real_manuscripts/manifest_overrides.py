"""Labels checked by hand (read by build_manifest.py).

Each entry was decided by reading lines of the transcription (and, where the repository
gives one, its description of the manuscript).  Where the repository's own table disagrees
with what the lines show, the note says so.
"""
OVERRIDE: dict[tuple[str, str], dict] = {}


def O(ds, ms, shelfmark="", form="", century="", language="", hand="", note=""):
    OVERRIDE[(ds, ms)] = {k: v for k, v in dict(shelfmark=shelfmark, form=form, century=century,
                                                language=language, hand=hand, note=note).items() if v}


G = "Gallicorpora-MSS-15"            # shelfmarks and titles from the TEI headers in data/*.xml
O(G, "btv1b100737746", "BnF NAF 6213", "prose", "15", note="Chronique fabuleuse (1180-1292)")
O(G, "btv1b10549431k", "BnF fr. 1376", "verse", "13-14", note="Florimont; Erec et Enide (TEI header: 1285-1315)")
O(G, "btv1b55008562q", "Arsenal 5103", "prose", "16", note="Bude, Institution du prince")
O(G, "btv1b60001038", "BnF fr. 12779", "verse", "15", note="Christine de Pizan")
O(G, "btv1b84260029", "BnF fr. 185", "prose", "", note="Legende des saints (Belet); no date in the TEI header")
O(G, "btv1b8448964p", "BnF fr. 604", "verse", "15", note="Christine de Pizan")
O(G, "btv1b90069505", "BnF fr. 1984", "prose", "15", note="Bertrand du Guesclin")
O(G, "btv1b90076543", "BnF fr. 413", "prose", "15", note="Vie des saints")
O(G, "btv1b9009670t", "BnF fr. 1419", "list", "15",
  note="Wauquelin, Alexandre: the five pages are the table of chapters (entries of 1-3 lines, each opening 'Coment')")
O(G, "btv1b90595325", "BnF fr. 777", "prose", "15", note="Berinus")
O(G, "btv1b90616967", "BnF fr. 12551", "prose", "15", note="Baudouin de Flandre")

C = "CREMMA-Medieval"                # centuries from the README table
O(C, "bnf_arsenal_3516-imageDuMonde", "Arsenal 3516", "verse", "13", note="Image du monde")
O(C, "bnf_fr_13496-saintJerome", "BnF fr. 13496", "prose", "13")
O(C, "bnf_fr_17229-saintLambert", "BnF fr. 17229", "prose", "13", note="other pages of this codex are in HTRomance-French")
O(C, "bnf_fr_1728", "BnF fr. 1728", "prose", "14")
O(C, "bnf_fr_22549-septSages", "BnF fr. 22549", "prose", "14", note="Sept sages de Rome")
O(C, "bnf_fr_24428-bestiaire", "BnF fr. 24428", "verse", "13", note="Bestiaire")
O(C, "bnf_fr_411-saintLambert", "BnF fr. 411", "prose", "14")
O(C, "bnf_fr_412-wauchier", "BnF fr. 412", "prose", "13", note="Wauchier de Denain, saints' lives")
O(C, "bnf_fr_844-manuscritDuRoi", "BnF fr. 844", "verse-run-on", "13", note="chansonnier with music; songs copied continuously")
O(C, "bodmer_168-otinel", "Cologny, Bodmer 168", "verse", "13", note="Otinel (laisses)")
O(C, "kbr_9232-examensMoraux", "KBR 9232", "prose", "15")
O(C, "pennsylvania_660-pelerinageMademoiselleSapience", "UPenn codex 660", "prose", "14")
O(C, "pennsylvania_codex_909-eneide", "UPenn codex 909", "verse", "15")
O(C, "vaticane_reg_lat_1616-otinel", "BAV Reg. lat. 1616", "verse", "14", note="Otinel (laisses)")

L = "CREMMA-Medieval-LAT"            # century from data-registry.csv; form from the lines
for ms in ("CCCC-MSS-236", "Latin6395", "Latin8236", "LaurentianusPluteus33.31", "LaurentianusPluteus39.34"):
    O(L, ms, form="verse")
for ms in ("Arras-861", "BGO-511", "BIS-193", "CCCC-MSS-165", "CLM13027", "H318", "Latin16195",
           "LaurentianusPluteus53.08", "LaurentianusPluteus53.09", "Mazarine915", "PalLat373", "Phi_10a135",
           "SBB_PK_Hdschr25", "UBL758", "WettF0015"):
    O(L, ms, form="prose")
O(L, "Egerton821", form="list", note="Sortes sanctorum: numbered lots, lines open with Roman numerals")

F = "HTRomance-French"
O(F, "bnf-fr-12581", form="verse-run-on", note="fabliau copied continuously (README: vers)")
O(F, "bnf-fr-1669", form="verse", note="README says prose; the lines are rhymed couplets")
O(F, "bnf-fr-12615", form="verse-run-on", note="chansonnier de Noailles")
O(F, "bnf-fr-20050", form="verse-run-on", note="chansonnier de Saint-Germain")
O(F, "bnf-fr-1881", form="verse", note="one hemistich per line")
O(F, "bnf-fr-13568", form="prose", century="14")
O(F, "bnf-fr-574", form="prose", century="14")
O(F, "bnf-fr-840", form="verse", century="14")
O(F, "bnf-arsenal-ms-3350", form="prose", century="15")

A = "HTRomance-Latin"
O(A, "bnf-lat-130", form="prose", note="README: mixed; Gregory, the lines are prose")
O(A, "bnf-lat-17903", form="verse", note="florilegium of verse extracts (README: mixed)")
O(A, "bnf-lat-17901", form="verse")
O(A, "bnf-arsenal-ms-1046", form="verse", note="Ovid")
O(A, "bnf-smith-lesouëf-11", form="verse")
for ms in ("bnf-lat-13388", "bnf-lat-15176", "bnf-lat-12449", "bnf-lat-12270", "bnf-nal-730", "bnf-lat-14650",
           "bnf-smith-lesouëf-12", "bnf-nal-632", "bnf-lat-6337", "bnf-smith-lesouëf-16"):
    O(A, ms, form="prose")

O("HTRomance-Italian", "bnf-it-912", form="list", note="price list of goods")
O("HTRomance-Spain", "bnf-esp-229", form="mixed", note="Laberinto de Fortuna with glosses")
O("HTRomance-Spain", "bnf-esp-98", language="Catalan")
O("HTRomance-Spain", "bnf-esp-548", language="Catalan")
O("HTRomance-Spain", "bnf-esp-225", language="Catalan")

P = "HTRomance-Occitan"
O(P, "paris-bnf-fr-14771", "BnF fr. 14771", "list", note="numbered lots (book of fortune)")
O(P, "paris-bnf-fr-19960", "BnF fr. 19960", "prose")
O(P, "paris-bnf-fr-25425", "BnF fr. 25425", "verse", "13", note="Canso de la Crozada (laisses); other pages of this codex are in HTRogene-Occitan")
O(P, "paris-bnf-fr-5728", "BnF fr. 5728", "prose")
O(P, "paris-bnf-fr-854", "BnF fr. 854", "verse-run-on", "13", note="chansonnier I: songs copied continuously, prose vidas")
O(P, "paris-bnf-naf-4232", form="exclude", note="duplicate: pages also in Gatmuza-Occitan bnf_naf_4232 (114 lines in common)")

HF = "HTRogene-French"
O(HF, "angers-am-cc03", form="documentary", hand="documentary", note="town accounts")
O(HF, "caen-ad14-g006", form="documentary", hand="documentary")
O(HF, "dijon-ad21-b11691", form="documentary", hand="documentary")
O(HF, "nantes-ad44-2mi31r1b1", form="documentary", hand="documentary", note="register")
O(HF, "paris-bnf-naf-4503", "BnF NAF 4503", "verse-run-on", "12", note="Vie de saint Alexis copied continuously (README: verse)")
O(HF, "paris-bnf-fr-146", "BnF fr. 146", "verse", "14", note="Roman de Fauvel")
O(HF, "paris-bnf-fr-12563", "BnF fr. 12563", "verse", "15")
O(HF, "paris-bnf-fr-12575", "BnF fr. 12575", "verse", "15")
O(HF, "paris-bnf-fr-345-346", form="list", note="one page, mostly a table of rubrics")
O(HF, "bordeaux-bm-731", form="prose")

HL = "HTRogene-Latin"
for ms in ("angers-archives-departementales-de-maine-et-loire-h-039-2-n-176-177",
           "angers-archives-departementales-de-maine-et-loire-h-045-1-ndeg-4",
           "auxerre-archives-departementales-de-l-yonne-h-2404", "besancon-bibliotheque-diocesaine-par-03",
           "laval-archives-de-la-mayenne-H-154", "liege-archives-de-l-etat-t51-12", "liege-archives-de-l-etat-t51-13",
           "liege-archives-de-l-etat-t51-14", "montigny-le-bretonneux-archives-departementales-des-yvelines-45h8-17",
           "montigny-le-bretonneux-archives-departementales-des-yvelines-45h8-8", "nice-am-aa-1-04",
           "paris-bibliotheque-de-l-ecole-nationale-superieure-des-beaux-arts-mn-mas-38",
           "pierrefitte-sur-seine-archives-nationales-de-france-ll-106b",
           "pierrefitte-sur-seine-archives-nationales-de-france-ll-108a",
           "pierrefitte-sur-seine-archives-nationales-de-france-ll-108b",
           "pierrefitte-sur-seine-archives-nationales-de-france-ll-110",
           "pierrefitte-sur-seine-archives-nationales-de-france-ll-125"):
    O(HL, ms, form="documentary", hand="documentary", note="charter, cartulary or chapter register")
O(HL, "paris-bibliotheque-nationale-de-france-lat-17226", "BnF lat. 17226", "prose", "7",
  note="uncial, written without word division: S3 not meaningful")
O(HL, "cambridge-corpus-christi-college-ms-29", "Cambridge, CCC 29", "prose", "13",
  note="Peter Comestor; first page is a table of chapters")

HI = "HTRogene-Italian"
O(HI, "oxford-bodleian-digby-94", form="prose", note="README says verse; the lines are prose")
O(HI, "strasbourg-bnu-1-808", form="prose", note="README says verse; the lines are prose")

HS = "HTRogene-Spanish"
O(HS, "paris-bnf-esp-43", "BnF esp. 43", "prose", "14", note="README says verse; long prose lines")
O(HS, "paris-bnf-esp-37", "BnF esp. 37", "verse", "15", note="cancionero")
O(HS, "paris-bnf-esp-305", "BnF esp. 305", "mixed", "15", language="Catalan")
O(HS, "paris-bnf-esp-479", form="verse")
O(HS, "rbme-y-iii-21", form="exclude", note="one usable page; other pages empty or uncorrected")
for ms in ("aca-general-n8", "barcelona-ahcb-deliberationes-1-b-ii-014"):
    O(HS, ms, form="documentary", hand="documentary", language="Catalan")
O(HS, "bc-ms-289", language="Catalan")
O(HS, "tarragona-ms-87-3", language="Catalan")

GO = "Gatmuza-Occitan"
O(GO, "Carcassonne_bm_34_flamenca", "Carcassonne BM 34", "verse", "14", note="Flamenca; other pages of this codex are in HTRogene-Occitan")
O(GO, "bnf_fr_12472_chansonnier_f", "BnF fr. 12472", "verse-run-on", "14")
O(GO, "bnf_fr_1592_chansonnier_b", "BnF fr. 1592", "verse-run-on", "13")
O(GO, "bnf_fr_2164_jaufre", "BnF fr. 2164", "verse", "13", note="Jaufre")
O(GO, "bnf_fr_2180_girart_de_rossilhon", "BnF fr. 2180", "verse", "13")
O(GO, "bnf_naf_4232_passion_didot_daurel", "BnF NAF 4232", "verse", "14")
O(GO, "montpellier_aa9_petit_thalamus", "Montpellier AM AA9", "documentary", "14", hand="documentary")

FB = "CIHAM-Fabliaux"
O(FB, "bern_354", "Bern, Burgerbibl. 354", "verse", "14")
O(FB, "bnf_Mazzarine_1553", form="exclude", note="duplicate: 281 of its 300 lines are the pages of HTRomance-French bnf-fr-1553 (BnF fr. 1553), transcribed a second time")
O(FB, "bnf_fr_1593", "BnF fr. 1593", "verse", "13")
O(FB, "bnf_fr_25545", "BnF fr. 25545", "verse", "14")
O(FB, "bnf_fr_837", "BnF fr. 837", "verse", "13")

LB = "CIHAM-Liber"
O(LB, "BnF_Naf_27401", "BnF NAF 27401", "prose", "14", note="Bersuire's Livy")
O(LB, "BnF_fr_263", "BnF fr. 263", "prose", "15", note="Bersuire's Livy")
O(LB, "Lisboa_bibl_nac_134", "Lisbon BN IL 134", "prose", "15", language="Latin", note="Trevet on Livy")
O(LB, "Lisboa_bibl_nac_135", "Lisbon BN IL 135", "prose", "15", language="Latin", note="Trevet on Livy")

O("DecameronFR", "DecameronFR", "Arsenal 5070", "prose", "15", note="Decameron in French (1430-1455)")
O("GrandesChroniques", "GrandesChroniques", "BnF fr. 2813", "prose", "14", note="Grandes Chroniques de France")
for t in "1234":
    O("Torino-L-II-14", f"IT-TO0265_MSL_II_14_Tome_0{t}", f"Turin BNU L.II.14, vol. {t}", "verse", "14",
      note="chansons de geste and romances; damaged by fire")
O("FoNDUE-LA-MSS-MA", "BnF_lat_7558", "BnF lat. 7558", "verse", "9", note="grammatical and poetic miscellany")
O("ParisBible", "ParisBible", "19 pages of 15 Paris Bibles (pooled)", "prose", "13",
  note="one or two pages per manuscript; pooled sample, not one manuscript")

O("Vienna-2023-Latin", "Vienna-2023-Latin", "ÖNB Cod. 3891", "prose", "15",
  note="Ebendorfer, Sermones (1441); abbreviations expanded by the transcribers")
O("Vienna-2024-Latin", "ONB_4135", "ÖNB Cod. 4135", "prose", "15", note="cursive; Tractatus contra Hussitas etc.")
O("Vienna-2024-Latin", "ONB_4680", "ÖNB Cod. 4680", "prose", "15", hand="autograph", note="Ebendorfer, autograph miscellany, cursive")
V5 = "Vienna-2025-Latin"
O(V5, "02_-_Praha,_NKP,_I_C_24_-_Winter_School", "Prague NKP I C 24", "prose", "14-15", note="Pulkava's Chronicle")
O(V5, "06_-_KNM,_sign__VIII_F_49_-_Winter_School", "Prague KNM VIII F 49", "prose", "14-15", note="Pulkava's Chronicle")
O(V5, "17_-_BSB,_Clm__476_-_Winter_School", "Munich BSB Clm 476", "prose", "14-15", note="Pulkava's Chronicle; includes annal lists ('Anno domini ...')")
O(V5, "19_-_Praha,_NKP,_sign__I_D_10_-_Winter_School", "Prague NKP I D 10", "prose", "14-15", note="Pulkava's Chronicle")
O(V5, "Aeneas_Silvius_Piccolomini_-_Wien,_Österreichische_Nationalbibliothek,_Cod__3445_-_75-196",
  "ÖNB Cod. 3445", "prose", "15", note="Aeneas Silvius; 25 of the 55 pages are marked 'under review' in the repository")
O(V5, "Chronicle_of_Dalimil_-_NKP_XII_E_17", "Prague NKP XII E 17", "verse", "14", note="Latin Dalimil fragment, one verse per line")
O(V5, "Cosmas_of_Prague_-_A3a_-_ÖNB_508_-_Winter_School", "ÖNB Cod. 508", "prose", "", note="Cosmas of Prague")

VG = "Vienna-2024-German"            # README names five ÖNB codices; file prefixes matched to them by content
O(VG, "DOD_51276", "ÖNB Cod. 3064 (matched by content)", "prose", "15", note="technical prose (gauging, gunnery)")
O(VG, "DTL_7978199", "ÖNB Cod. 2875 (matched by content)", "prose", "15", note="theological miscellany c. 1400, narrow columns")
O(VG, "DTL_9133478", "ÖNB Cod. 2780 (matched by content)", "list", "15",
  note="law book (1423): the transcribed pages are lists of chapter titles, one title per line")
O(VG, "Renner Cod. 2810", "ÖNB Cod. 2810", "verse", "15", note="Hugo von Trimberg, Der Renner (1410)")
O(VG, "enikel weltchronik", "ÖNB Cod. 2921", "verse", "14", note="Jans Enikel, Weltchronik (1397)")

O("Vienna-2023-Czech", "Vienna-2023-Czech", "ÖNB Cod. 3304", "prose", "15",
  note="New Testament of Martin Lupac (1440); abbreviations expanded")
O("Vienna-2025-Czech", "Vienna-2025-Czech", "Krakow BJ Rkp. 441 IV", "prose", "16", note="Czech chronicle (1531)")
O("CREMMA-MSS-16", "CREMMA-MSS-16", "9 pages of 16th-century French documents (pooled)", "documentary", "16", hand="documentary")
O("Vienna-2025-Hebrew", "Vienna-2025-Hebrew", "BL Or. 10037", "prose", "15", note="Mishneh Torah; right-to-left, geometry mirrored")
O("Vienna-2024-Syriac", "Vienna-2024-Syriac", "ÖNB Cod. Syr. 1", "prose", "16", note="right-to-left, geometry mirrored")
O("Vienna-2025-Syriac", "Vienna-2025-Syriac", "Jerusalem SMMJ 36", "prose", "12-14", note="right-to-left, geometry mirrored")
O("LiDi-Ligorio", "LiDi-Ligorio", "Pirro Ligorio autograph (files Ja_III)", "list", "16", hand="autograph",
  note="antiquarian notes with copied inscriptions; one block per entry")
O("CarolineMinuscule", "CarolineMinuscule", "17 pages of 14 BSB manuscripts (pooled)", "prose", "9-12",
  note="pooled sample, not one manuscript")
O("iForal", "iForal", "Portuguese forais (pooled)", "documentary", "13-15", hand="documentary",
  note="municipal charters, pooled")
O("GasparoSardi", "GasparoSardi", "Bern, Burgerbibl. Cod. 174", "list", "16", hand="autograph",
  note="Gasparo Sardi, Toponomasia: alphabetical dictionary of place names")
O("BibliaArabica", "OxfordLaudOr258", "Oxford, Bodl. Laud Or. 258", "prose", "", note="right-to-left, geometry mirrored")
