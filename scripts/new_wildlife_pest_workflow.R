# =============================================================================
# BA SPECIES WILDLIFE TRAITS & PEST/PATHOGEN WORKFLOW
# Trees only — Lower 48 United States
# Zero hardcoded trait values — all data sourced from downloadable databases
# =============================================================================
#
# DATA SOURCES — what to download and where:
#
#  1. plants.db          — already have (plantatlas.ai)
#  2. armycorps_wetland_plants.xlsx — already have (Army Corps)
#  3. HOSTS_lep.csv      — Natural History Museum London, CC0 public domain
#                          Download: https://data.nhm.ac.uk/dataset/hosts
#                          Click "HOSTS" resource -> Download (140,485 records)
#                          Save as: HOSTS_lep.csv
#  4. NIDRM_pests.xlsx   — USDA Forest Service, public domain
#                          Download: https://www.fs.usda.gov/foresthealth/
#                            applied-sciences/mapping-reporting/national-risk-maps
#                          Click "2013-2027 NIDRM Model Details workbook"
#                          Save as: NIDRM_pests.xlsx
#
# WHAT IS STILL HARDCODED AND WHY IT CANNOT BE SOURCED ELSEWHERE:
#
#  - sci_name_map (BA code -> scientific name): this is a BA-specific
#    translation table. No external database knows what "quercusRubra" means.
#    It MUST be hardcoded. It is not trait data — it is a key lookup.
#
#  - NON_SPECIES list: same reason — BA-specific metadata.
#
# EVERYTHING ELSE comes from the four downloaded files above.
# =============================================================================

library(tidyverse)
library(janitor)
library(RSQLite)
library(DBI)
library(readxl)

# =============================================================================
# STEP 0: FILE PATHS
# =============================================================================
PATH_PLANTS_DB  <- "~/Documents/PBGJAM-data-explorer/raw_data/v2/plants.db"
PATH_WETLAND    <- "~/Documents/PBGJAM-data-explorer/raw_data/v2/armycorps_wetland_plants.xlsx"
PATH_HOSTS_LEP <- "~/Documents/PBGJAM-data-explorer/raw_data/v2/HOSTS_lep.csv"       # NHM London Lepidoptera hostplant DB
PATH_NIDRM      <- "~/Documents/PBGJAM-data-explorer/raw_data/v2/NIDRM_pests.xlsx"    # USFS National Insect & Disease Risk Map workbook
PATH_BA         <- "~/Documents/PBGJAM-data-explorer/data4Tate/BA.rdata"
# =============================================================================
# STEP 1: LOAD BA AND BUILD SPECIES LOOKUP
# (sci_name_map is the ONLY hardcoded content — it is a key, not trait data)
# =============================================================================
load(PATH_BA)
raw_species <- colnames(BA)

NON_SPECIES <- c("family ANot list", "treeBroadlea", "treeEvergree", "treeUnknown")

parse_ba_name <- function(x) {
  tibble(
    ba_code   = x,
    genus     = str_to_title(str_extract(x, "^[a-z]+")),
    sp_abbrev = str_extract(x, "(?<=[a-z])[A-Z][a-zA-Z]+$")
  )
}

sci_name_map <- tribble(
  ~ba_code,               ~scientific_name,
  "abiesAmabilis",        "Abies amabilis",
  "abiesBalsamea",        "Abies balsamea",
  "abiesConcolor",        "Abies concolor",
  "abiesFraseri",         "Abies fraseri",
  "abiesGrandis",         "Abies grandis",
  "abiesLasiocar",        "Abies lasiocarpa",
  "abiesMagnific",        "Abies magnifica",
  "abiesProcera",         "Abies procera",
  "abiesShastens",        "Abies shastensis",
  "acaciaConfusa",        "Acacia confusa",
  "acaciaFarnesia",       "Acacia farnesiana",
  "acaciaKoa",            "Acacia koa",
  "acaciaKoaia",          "Acacia koaia",
  "acerBarbatum",         "Acer barbatum",
  "acerGlabrum",          "Acer glabrum",
  "acerGrandide",         "Acer grandidentatum",
  "acerLeucoder",         "Acer leucoderme",
  "acerMacrophy",         "Acer macrophyllum",
  "acerNegundo",          "Acer negundo",
  "acerNigrum",           "Acer nigrum",
  "acerPensylva",         "Acer pensylvanicum",
  "acerPlatanoi",         "Acer platanoides",
  "acerRubrum",           "Acer rubrum",
  "acerSacchari",         "Acer saccharinum",
  "acerSaccharu",         "Acer saccharum",
  "acerSpicatum",         "Acer spicatum",
  "aesculusCaliforn",     "Aesculus californica",
  "aesculusFlava",        "Aesculus flava",
  "aesculusGlabra",       "Aesculus glabra",
  "ailanthuAltissim",     "Ailanthus altissima",
  "albiziaChinensi",      "Albizia chinensis",
  "albiziaJulibris",      "Albizia julibrissin",
  "aleuriteMoluccan",     "Aleurites moluccanus",
  "alnusGlutinos",        "Alnus glutinosa",
  "alnusOblongif",        "Alnus oblongifolia",
  "alnusRhombifo",        "Alnus rhombifolia",
  "alnusRubra",           "Alnus rubra",
  "amelanchArborea",      "Amelanchier arborea",
  "annonaGlabra",         "Annona glabra",
  "arbutusArizonic",      "Arbutus arizonica",
  "arbutusMenziesi",      "Arbutus menziesii",
  "arbutusXalapens",      "Arbutus xalapensis",
  "asiminaTriloba",       "Asimina triloba",
  "avicenniGerminan",     "Avicennia germinans",
  "betulaAlleghan",       "Betula alleghaniensis",
  "betulaLenta",          "Betula lenta",
  "betulaNeoalask",       "Betula neoalaskana",
  "betulaNigra",          "Betula nigra",
  "betulaOccident",       "Betula occidentalis",
  "betulaPapyrife",       "Betula papyrifera",
  "betulaPopulifo",       "Betula populifolia",
  "burseraSimaruba",      "Bursera simaruba",
  "calocedrDecurren",     "Calocedrus decurrens",
  "carpinusCarolini",     "Carpinus caroliniana",
  "caryaAlba",            "Carya alba",
  "caryaAquatica",        "Carya aquatica",
  "caryaCarolina",        "Carya carolinae-septentrionalis",
  "caryaCordifor",        "Carya cordiformis",
  "caryaFloridan",        "Carya floridana",
  "caryaGlabra",          "Carya glabra",
  "caryaIllinoin",        "Carya illinoinensis",
  "caryaLacinios",        "Carya laciniosa",
  "caryaMyristic",        "Carya myristiciformis",
  "caryaOvalis",          "Carya ovalis",
  "caryaOvata",           "Carya ovata",
  "caryaPallida",         "Carya pallida",
  "caryaTexana",          "Carya texana",
  "castaneaDentata",      "Castanea dentata",
  "castaneaMollissi",     "Castanea mollissima",
  "castaneaPumila",       "Castanea pumila",
  "casuarinCunningh",     "Casuarina cunninghamiana",
  "casuarinEquiseti",     "Casuarina equisetifolia",
  "casuarinGlauca",       "Casuarina glauca",
  "casuarinLepidoph",     "Casuarina lepidophloia",
  "catalpaBignonio",      "Catalpa bignonioides",
  "catalpaSpeciosa",      "Catalpa speciosa",
  "celtisLaevigat",       "Celtis laevigata",
  "celtisOccident",       "Celtis occidentalis",
  "cercisCanadens",       "Cercis canadensis",
  "cercocarLedifoli",     "Cercocarpus ledifolius",
  "chamaecyLawsonia",     "Chamaecyparis lawsoniana",
  "chamaecyNootkate",     "Chamaecyparis nootkatensis",
  "chamaecyThyoides",     "Chamaecyparis thyoides",
  "cinnamomCamphora",     "Cinnamomum camphora",
  "cladrastKentukea",     "Cladrastis kentukea",
  "cocosNucifera",        "Cocos nucifera",
  "condaliaHookeri",      "Condalia hookeri",
  "conocarpErectus",      "Conocarpus erectus",
  "cornusFlorida",        "Cornus florida",
  "cornusNuttalli",       "Cornus nuttallii",
  "cotinusObovatus",      "Cotinus obovatus",
  "crataeguCrus-gal",     "Crataegus crus-galli",
  "crataeguMollis",       "Crataegus mollis",
  "cupressuArizonic",     "Cupressus arizonica",
  "cupressuBakeri",       "Cupressus bakeri",
  "cupressuMacrocar",     "Cupressus macrocarpa",
  "cupressuSargenti",     "Cupressus sargentii",
  "diospyroTexana",       "Diospyros texana",
  "diospyroVirginia",     "Diospyros virginiana",
  "ebenopsiEbano",        "Ebenopsis ebano",
  "ehretiaAnacua",        "Ehretia anacua",
  "elaeagnuAngustif",     "Elaeagnus angustifolia",
  "eucalyptGlobulus",     "Eucalyptus globulus",
  "eucalyptGrandis",      "Eucalyptus grandis",
  "eucalyptRobusta",      "Eucalyptus robusta",
  "eucalyptSaligna",      "Eucalyptus saligna",
  "fagusGrandifo",        "Fagus grandifolia",
  "ficusAurea",           "Ficus aurea",
  "ficusMicrocar",        "Ficus microcarpa",
  "fraxinusAmerican",     "Fraxinus americana",
  "fraxinusCarolini",     "Fraxinus caroliniana",
  "fraxinusLatifoli",     "Fraxinus latifolia",
  "fraxinusNigra",        "Fraxinus nigra",
  "fraxinusPennsylv",     "Fraxinus pennsylvanica",
  "fraxinusProfunda",     "Fraxinus profunda",
  "fraxinusQuadrang",     "Fraxinus quadrangulata",
  "fraxinusTexensis",     "Fraxinus texensis",
  "fraxinusUhdei",        "Fraxinus uhdei",
  "fraxinusVelutina",     "Fraxinus velutina",
  "ginkgoBiloba",         "Ginkgo biloba",
  "gleditsiAquatica",     "Gleditsia aquatica",
  "gleditsiTriacant",     "Gleditsia triacanthos",
  "gordoniaLasianth",     "Gordonia lasianthus",
  "grevilleRobusta",      "Grevillea robusta",
  "gymnoclaDioicus",      "Gymnocladus dioicus",
  "halesiaCarolina",      "Halesia carolina",
  "halesiaDiptera",       "Halesia diptera",
  "hibiscusTiliaceu",     "Hibiscus tiliaceus",
  "ilexOpaca",            "Ilex opaca",
  "juglansCaliforn",      "Juglans californica",
  "juglansCinerea",       "Juglans cinerea",
  "juglansHindsii",       "Juglans hindsii",
  "juglansMajor",         "Juglans major",
  "juglansMicrocar",      "Juglans microcarpa",
  "juglansNigra",         "Juglans nigra",
  "juniperuAshei",        "Juniperus ashei",
  "juniperuCaliforn",     "Juniperus californica",
  "juniperuCoahuile",     "Juniperus coahuilensis",
  "juniperuDeppeana",     "Juniperus deppeana",
  "juniperuFlaccida",     "Juniperus flaccida",
  "juniperuMonosper",     "Juniperus monosperma",
  "juniperuOccident",     "Juniperus occidentalis",
  "juniperuOsteospe",     "Juniperus osteosperma",
  "juniperuPinchoti",     "Juniperus pinchotii",
  "juniperuScopulor",     "Juniperus scopulorum",
  "juniperuVirginia",     "Juniperus virginiana",
  "lagunculRacemosa",     "Laguncularia racemosa",
  "larixLaricina",        "Larix laricina",
  "larixLyallii",         "Larix lyallii",
  "larixOccident",        "Larix occidentalis",
  "leucaenaLeucocep",     "Leucaena leucocephala",
  "liquidamStyracif",     "Liquidambar styraciflua",
  "liriodenTulipife",     "Liriodendron tulipifera",
  "lithocarDensiflo",     "Lithocarpus densiflorus",
  "macluraPomifera",      "Maclura pomifera",
  "magnoliaAcuminat",     "Magnolia acuminata",
  "magnoliaFraseri",      "Magnolia fraseri",
  "magnoliaGrandifl",     "Magnolia grandiflora",
  "magnoliaMacrophy",     "Magnolia macrophylla",
  "magnoliaTripetal",     "Magnolia tripetala",
  "magnoliaVirginia",     "Magnolia virginiana",
  "malusAngustif",        "Malus angustifolia",
  "malusCoronari",        "Malus coronaria",
  "malusFusca",           "Malus fusca",
  "malusIoensis",         "Malus ioensis",
  "meliaAzedarac",        "Melia azedarach",
  "metopiumToxiferu",     "Metopium toxiferum",
  "morellaFaya",          "Morella faya",
  "morusAlba",            "Morus alba",
  "morusMicrophy",        "Morus microphylla",
  "morusRubra",           "Morus rubra",
  "nyssaAquatica",        "Nyssa aquatica",
  "nyssaBiflora",         "Nyssa biflora",
  "nyssaOgeche",          "Nyssa ogeche",
  "nyssaSylvatic",        "Nyssa sylvatica",
  "olneyaTesota",         "Olneya tesota",
  "ostryaVirginia",       "Ostrya virginiana",
  "oxydendrArboreum",     "Oxydendrum arboreum",
  "paulowniTomentos",     "Paulownia tomentosa",
  "perseaAmerican",       "Persea americana",
  "perseaBorbonia",       "Persea borbonia",
  "piceaAbies",           "Picea abies",
  "piceaBreweria",        "Picea breweriana",
  "piceaEngelman",        "Picea engelmannii",
  "piceaGlauca",          "Picea glauca",
  "piceaMariana",         "Picea mariana",
  "piceaPungens",         "Picea pungens",
  "piceaRubens",          "Picea rubens",
  "piceaSitchens",        "Picea sitchensis",
  "pinusAlbicaul",        "Pinus albicaulis",
  "pinusAristata",        "Pinus aristata",
  "pinusArizonic",        "Pinus arizonica",
  "pinusAttenuat",        "Pinus attenuata",
  "pinusBalfouri",        "Pinus balfouriana",
  "pinusBanksian",        "Pinus banksiana",
  "pinusCembroid",        "Pinus cembroides",
  "pinusClausa",          "Pinus clausa",
  "pinusContorta",        "Pinus contorta",
  "pinusCoulteri",        "Pinus coulteri",
  "pinusDiscolor",        "Pinus discolor",
  "pinusEchinata",        "Pinus echinata",
  "pinusEdulis",          "Pinus edulis",
  "pinusElliotti",        "Pinus elliottii",
  "pinusEngelman",        "Pinus engelmannii",
  "pinusFlexilis",        "Pinus flexilis",
  "pinusGlabra",          "Pinus glabra",
  "pinusJeffreyi",        "Pinus jeffreyi",
  "pinusLamberti",        "Pinus lambertiana",
  "pinusLeiophyl",        "Pinus leiophylla",
  "pinusLongaeva",        "Pinus longaeva",
  "pinusMonophyl",        "Pinus monophylla",
  "pinusMonticol",        "Pinus monticola",
  "pinusMuricata",        "Pinus muricata",
  "pinusNigra",           "Pinus nigra",
  "pinusPalustri",        "Pinus palustris",
  "pinusPatula",          "Pinus patula",
  "pinusPinaster",        "Pinus pinaster",
  "pinusPonderos",        "Pinus ponderosa",
  "pinusPungens",         "Pinus pungens",
  "pinusRadiata",         "Pinus radiata",
  "pinusRemota",          "Pinus remota",
  "pinusResinosa",        "Pinus resinosa",
  "pinusRigida",          "Pinus rigida",
  "pinusSabinian",        "Pinus sabiniana",
  "pinusSerotina",        "Pinus serotina",
  "pinusStrobifo",        "Pinus strobiformis",
  "pinusStrobus",         "Pinus strobus",
  "pinusSylvestr",        "Pinus sylvestris",
  "pinusTaeda",           "Pinus taeda",
  "pinusVirginia",        "Pinus virginiana",
  "pinusWashoens",        "Pinus washoensis",
  "piscidiaPiscipul",     "Piscidia piscipula",
  "platanusOccident",     "Platanus occidentalis",
  "platanusRacemosa",     "Platanus racemosa",
  "platanusWrightii",     "Platanus wrightii",
  "populusAlba",          "Populus alba",
  "populusAngustif",      "Populus angustifolia",
  "populusBalsamif",      "Populus balsamifera",
  "populusDeltoide",      "Populus deltoides",
  "populusFremonti",      "Populus fremontii",
  "populusGrandide",      "Populus grandidentata",
  "populusHeteroph",      "Populus heterophylla",
  "populusTremuloi",      "Populus tremuloides",
  "prosopisGlandulo",     "Prosopis glandulosa",
  "prosopisPallida",      "Prosopis pallida",
  "prosopisPubescen",     "Prosopis pubescens",
  "prosopisVelutina",     "Prosopis velutina",
  "prunusAmerican",       "Prunus americana",
  "prunusAvium",          "Prunus avium",
  "prunusCerasus",        "Prunus cerasus",
  "prunusEmargina",       "Prunus emarginata",
  "prunusNigra",          "Prunus nigra",
  "prunusPensylva",       "Prunus pensylvanica",
  "prunusPersica",        "Prunus persica",
  "prunusSerotina",       "Prunus serotina",
  "prunusVirginia",       "Prunus virginiana",
  "pseudotsMacrocar",     "Pseudotsuga macrocarpa",
  "pseudotsMenziesi",     "Pseudotsuga menziesii",
  "psidiumCattleia",      "Psidium cattleianum",
  "psidiumGuajava",       "Psidium guajava",
  "quercusAgrifoli",      "Quercus agrifolia",
  "quercusAlba",          "Quercus alba",
  "quercusArizonic",      "Quercus arizonica",
  "quercusBicolor",       "Quercus bicolor",
  "quercusBuckleyi",      "Quercus buckleyi",
  "quercusChrysole",      "Quercus chrysolepis",
  "quercusCoccinea",      "Quercus coccinea",
  "quercusDouglasi",      "Quercus douglasii",
  "quercusEllipsoi",      "Quercus ellipsoidalis",
  "quercusEmoryi",        "Quercus emoryi",
  "quercusEngelman",      "Quercus engelmannii",
  "quercusFalcata",       "Quercus falcata",
  "quercusGambelii",      "Quercus gambelii",
  "quercusGarryana",      "Quercus garryana",
  "quercusGracilif",      "Quercus graciliformis",
  "quercusGravesii",      "Quercus gravesii",
  "quercusGrisea",        "Quercus grisea",
  "quercusHypoleuc",      "Quercus hypoleucoides",
  "quercusIlicifol",      "Quercus ilicifolia",
  "quercusImbricar",      "Quercus imbricaria",
  "quercusIncana",        "Quercus incana",
  "quercusKelloggi",      "Quercus kelloggii",
  "quercusLaceyi",        "Quercus laceyi",
  "quercusLaevis",        "Quercus laevis",
  "quercusLaurifol",      "Quercus laurifolia",
  "quercusLobata",        "Quercus lobata",
  "quercusLyrata",        "Quercus lyrata",
  "quercusMacrocar",      "Quercus macrocarpa",
  "quercusMargaret",      "Quercus margaretta",
  "quercusMariland",      "Quercus marilandica",
  "quercusMichauxi",      "Quercus michauxii",
  "quercusMinima",        "Quercus minima",
  "quercusMuehlenb",      "Quercus muehlenbergii",
  "quercusNigra",         "Quercus nigra",
  "quercusOblongif",      "Quercus oblongifolia",
  "quercusOglethor",      "Quercus oglethorpensis",
  "quercusPagoda",        "Quercus pagoda",
  "quercusPalustri",      "Quercus palustris",
  "quercusPhellos",       "Quercus phellos",
  "quercusPrinoide",      "Quercus prinoides",
  "quercusPrinus",        "Quercus prinus",
  "quercusRubra",         "Quercus rubra",
  "quercusRugosa",        "Quercus rugosa",
  "quercusShumardi",      "Quercus shumardii",
  "quercusSimilis",       "Quercus similis",
  "quercusSinuata",       "Quercus sinuata",
  "quercusStellata",      "Quercus stellata",
  "quercusTexana",        "Quercus texana",
  "quercusVelutina",      "Quercus velutina",
  "quercusVirginia",      "Quercus virginiana",
  "quercusWislizen",      "Quercus wislizeni",
  "rhizophoMangle",       "Rhizophora mangle",
  "robiniaNeomexic",      "Robinia neomexicana",
  "robiniaPseudoac",      "Robinia pseudoacacia",
  "sabalPalmetto",        "Sabal palmetto",
  "salixAlba",            "Salix alba",
  "salixAmygdalo",        "Salix amygdaloides",
  "salixBebbiana",        "Salix bebbiana",
  "salixCarolini",        "Salix caroliniana",
  "salixNigra",           "Salix nigra",
  "salixPyrifoli",        "Salix pyrifolia",
  "salixX sepulc",        "Salix x sepulcralis",
  "samaneaSaman",         "Samanea saman",
  "sapindusSaponari",     "Sapindus saponaria",
  "sassafraAlbidum",      "Sassafras albidum",
  "schinusTerebint",      "Schinus terebinthifolia",
  "sequoiadGiganteu",     "Sequoiadendron giganteum",
  "sequoiaSempervi",      "Sequoia sempervirens",
  "sideroxyLanugino",     "Sideroxylon lanuginosum",
  "sideroxySalicifo",     "Sideroxylon salicifolium",
  "sophoraAffinis",       "Sophora affinis",
  "sorbusAmerican",       "Sorbus americana",
  "sorbusAucupari",       "Sorbus aucuparia",
  "sorbusDecora",         "Sorbus decora",
  "spathodeCampanul",     "Spathodea campanulata",
  "swieteniMacrophy",     "Swietenia macrophylla",
  "swieteniMahagoni",     "Swietenia mahagoni",
  "syzygiumCumini",       "Syzygium cumini",
  "syzygiumJambos",       "Syzygium jambos",
  "tamarixUNKN",          "Tamarix sp.",
  "taxodiumAscenden",     "Taxodium ascendens",
  "taxodiumDistichu",     "Taxodium distichum",
  "taxusBrevifol",        "Taxus brevifolia",
  "thujaOccident",        "Thuja occidentalis",
  "thujaPlicata",         "Thuja plicata",
  "tiliaAmerican",        "Tilia americana",
  "toonaCiliata",         "Toona ciliata",
  "torreyaCaliforn",      "Torreya californica",
  "tremaOriental",        "Trema orientale",
  "triadicaSebifera",     "Triadica sebifera",
  "tsugaCanadens",        "Tsuga canadensis",
  "tsugaCarolini",        "Tsuga caroliniana",
  "tsugaHeteroph",        "Tsuga heterophylla",
  "tsugaMertensi",        "Tsuga mertensiana",
  "ulmusAlata",           "Ulmus alata",
  "ulmusAmerican",        "Ulmus americana",
  "ulmusCrassifo",        "Ulmus crassifolia",
  "ulmusPumila",          "Ulmus pumila",
  "ulmusRubra",           "Ulmus rubra",
  "ulmusSerotina",        "Ulmus serotina",
  "ulmusThomasii",        "Ulmus thomasii",
  "umbellulCaliforn",     "Umbellularia californica",
  "verniciaFordii",       "Vernicia fordii",
  # ── Corrected names based on plants.db lookup ──────────────────────────────
  # Reclassifications and synonyms confirmed in DB
  "acaciaFarnesia",       "Vachellia farnesiana",       # Acacia -> Vachellia
  "acerBarbatum",         "Acer floridanum",            # barbatum is synonym
  "caryaAlba",            "Carya tomentosa",            # alba is synonym of tomentosa
  "lithocarDensiflo",     "Notholithocarpus densiflorus", # reclassified genus
  "quercusPrinus",        "Quercus montana",            # prinus is synonym of montana
  "schinusTerebint",      "Schinus terebinthifolius",   # DB spelling: folius not folia
  "tremaOriental",        "Trema micrantha",            # orientale -> micrantha in DB
  # Cupressus spp. are in DB as Hesperocyparis (USDA reclassification)
  "cupressuArizonic",     "Hesperocyparis arizonica",
  "cupressuBakeri",       "Hesperocyparis bakeri",
  "cupressuMacrocar",     "Hesperocyparis macrocarpa",
  "cupressuSargenti",     "Hesperocyparis sargentii",
  "cupressuMacnabia",     "Hesperocyparis macnabiana",
  # Chamaecyparis nootkatensis -> in DB as Callitropsis nootkatensis
  "chamaecyNootkate",     "Callitropsis nootkatensis",
  # Sophora affinis -> reclassified as Styphnolobium affine
  "sophoraAffinis",       "Styphnolobium affine",
  "chrysoleChrysoph",     "Chrysolepis chrysophylla"  # golden chinkapin - CA native tree
)

species_df <- tibble(ba_code = raw_species) %>%
  filter(!ba_code %in% NON_SPECIES) %>%
  mutate(is_unknown = str_detect(ba_code, "UNKN")) %>%
  mutate(parse_ba_name(ba_code)) %>%
  left_join(sci_name_map, by = "ba_code") %>%
  mutate(
    scientific_name = if_else(is.na(scientific_name),
                              paste(genus, str_to_lower(sp_abbrev)),
                              scientific_name),
    name_verified = ba_code %in% sci_name_map$ba_code
  )

message(sprintf("Species parsed: %d total | %d unknowns | %d verified names",
                nrow(species_df), sum(species_df$is_unknown),
                sum(species_df$name_verified)))


# =============================================================================
# STEP 2: QUERY plants.db — TREES ONLY
# Filter to growth_habits containing "Tree" to exclude shrubs and grasses.
# growth_habits lives in the `plants` table (e.g. "Tree", "Tree, Shrub").
# =============================================================================
con <- dbConnect(SQLite(), PATH_PLANTS_DB)

# scientific_name in this DB is HTML-wrapped with author strings, e.g.:
#   "<i>Quercus rubra</i> L."
# We must clean in R (SQLite has no regex replace), then match.
# Also: each binomial appears as Species + Varieties/Subspecies rows.
# We keep only rank = "Species" to get one row per binomial.

all_plants_raw <- dbGetQuery(con,
                             "SELECT id, symbol, scientific_name, common_name,
          growth_habits, durations, rank
   FROM plants
   WHERE rank = 'Species'
     AND (growth_habits LIKE '%Tree%' OR growth_habits IS NULL)"
) %>% as_tibble() %>%
  mutate(
    clean_scientific_name = scientific_name %>%
      str_remove_all("<[^>]+>") %>%   # strip <i>...</i> and any other tags
      str_trim() %>%
      str_extract("^\\S+\\s+\\S+") %>%  # keep only Genus + epithet
      str_trim()
  )

message(sprintf("plants table (Species rank, tree habit): %d rows, %d unique clean names",
                nrow(all_plants_raw),
                n_distinct(all_plants_raw$clean_scientific_name)))

# Join to our species list on clean name
id_lookup <- species_df %>%
  select(scientific_name) %>%
  inner_join(
    all_plants_raw %>%
      select(plant_id = id, symbol, clean_scientific_name,
             common_name, growth_habits, durations),
    by = c("scientific_name" = "clean_scientific_name")
  ) %>%
  distinct(scientific_name, .keep_all = TRUE)

message(sprintf("plants table: %d of %d species matched", nrow(id_lookup), nrow(species_df)))

if (nrow(id_lookup) == 0) {
  stop("No species matched. Check that species_df$scientific_name uses plain binomials (Genus species).")
}

# Pull characteristics for matched plant_ids
if (nrow(id_lookup) > 0) {
  plant_ids_sql <- paste(id_lookup$plant_id, collapse = ",")
  
  chars_raw <- dbGetQuery(con, sprintf(
    "SELECT
       plant_id,
       active_growth_period,
       growth_rate,
       growth_form,
       height_mature_feet,
       lifespan,
       resprout_ability,
       seedling_vigor,
       shade_tolerance,
       drought_tolerance,
       moisture_use,
       anaerobic_tolerance,
       fire_tolerance,
       adapted_to_coarse_textured_soils,
       adapted_to_fine_textured_soils,
       adapted_to_medium_textured_soils,
       ph_minimum,
       ph_maximum,
       salinity_tolerance,
       caco_3_tolerance,
       berry_nut_seed_product,
       fruit_seed_abundance,
       fruit_seed_period_begin,
       fruit_seed_period_end,
       fruit_seed_persistence,
       palatable_browse_animal,
       palatable_human,
       bloom_period,
       nitrogen_fixation,
       known_allelopath,
       seed_spread_rate,
       vegetative_spread_rate,
       toxicity
     FROM plant_characteristics
     WHERE plant_id IN (%s)",
    plant_ids_sql
  )) %>% as_tibble() %>% clean_names()
  
  message(sprintf("plant_characteristics: %d rows returned", nrow(chars_raw)))
  
  # growth_form is sparsely populated — most tree records have NULL.
  # The main tree filter already happened via growth_habits in the plants query.
  # Here we only exclude rows explicitly flagged as non-tree growth forms.
  # Keeping NULLs is correct: most trees in this DB have no growth_form entry.
  non_tree_forms <- c("shrub", "subshrub", "graminoid", "forb",
                      "vine", "nonvascular", "lichenous")
  
  tree_form_chars <- chars_raw %>%
    filter(
      is.na(growth_form) |
        !str_detect(tolower(growth_form), paste(non_tree_forms, collapse = "|"))
    )
  
  message(sprintf("plant_characteristics: %d rows kept | %d excluded as non-tree",
                  nrow(tree_form_chars), nrow(chars_raw) - nrow(tree_form_chars)))
  
  if (nrow(tree_form_chars) == 0) {
    warning("tree_form_chars is empty after filter. growth_form values in chars_raw: ",
            paste(unique(chars_raw$growth_form), collapse = ", "))
  }
  
  # Native status
  native_raw <- dbGetQuery(con, sprintf(
    "SELECT plant_id, region, status, status_type
     FROM plant_native_statuses
     WHERE plant_id IN (%s)",
    plant_ids_sql
  )) %>% as_tibble()
  
  native_summary <- native_raw %>%
    group_by(plant_id) %>%
    summarise(
      native_l48      = any(region == "L48" & status_type == "Native"),
      introduced_l48  = any(region == "L48" & status_type == "Introduced"),
      native_regions  = paste(region[status_type == "Native"], collapse = ","),
      .groups = "drop"
    )
  
  distrib_raw <- dbGetQuery(con, sprintf(
    "SELECT plant_id, native_states, introduced_states
     FROM plant_distribution WHERE plant_id IN (%s)",
    plant_ids_sql
  )) %>% as_tibble()
  
} else {
  tree_form_chars <- native_summary <- distrib_raw <- tibble()
}

dbDisconnect(con)

# Use left_join so species in id_lookup without characteristics data are kept
# (rather than silently dropped by inner_join when tree_form_chars is sparse)
plants_clean <- id_lookup %>%
  left_join(tree_form_chars, by = "plant_id") %>%
  left_join(native_summary,  by = "plant_id") %>%
  left_join(distrib_raw,     by = "plant_id") %>%
  rename(soil_ph_min = ph_minimum, soil_ph_max = ph_maximum,
         height_mature_ft = height_mature_feet) %>%
  distinct(scientific_name, .keep_all = TRUE)

message(sprintf("plants_clean: %d species | %d with characteristics data",
                nrow(plants_clean),
                sum(!is.na(plants_clean$growth_rate))))

message(sprintf("Final tree species with plant data: %d", nrow(plants_clean)))

# =============================================================================
# STEP 3: ARMY CORPS WETLAND — LOWER 48
# =============================================================================
wetland_raw <- read_excel(PATH_WETLAND, sheet = 1) %>%
  clean_names() %>%
  rename_with(~ "scientific_name", 1) %>%
  mutate(scientific_name = str_trim(scientific_name))

L48_COLS <- c("agcp", "aw", "emp", "gp", "mw", "ncne", "wmvc")
avail_l48 <- intersect(L48_COLS, names(wetland_raw))

wetland_rank <- c("OBL" = 5, "FACW" = 4, "FAC" = 3, "FACU" = 2, "UPL" = 1)

wetland_l48 <- wetland_raw %>%
  select(scientific_name, all_of(avail_l48)) %>%
  mutate(across(all_of(avail_l48), ~ na_if(trimws(.x), ""))) %>%
  distinct(scientific_name, .keep_all = TRUE) %>%
  mutate(
    wetland_wettest = apply(select(., all_of(avail_l48)), 1, function(r) {
      v <- na.omit(r); if (!length(v)) return(NA_character_)
      v[which.max(wetland_rank[v])]
    }),
    wetland_modal = apply(select(., all_of(avail_l48)), 1, function(r) {
      v <- na.omit(r); if (!length(v)) return(NA_character_)
      names(sort(table(v), decreasing = TRUE))[1]
    }),
    riparian_category = case_when(
      wetland_wettest == "OBL"  ~ "Obligate Wetland",
      wetland_wettest == "FACW" ~ "Facultative Wetland — riparian preferred",
      wetland_wettest == "FAC"  ~ "Facultative — tolerates riparian",
      wetland_wettest == "FACU" ~ "Facultative Upland — rarely riparian",
      wetland_wettest == "UPL"  ~ "Upland Only",
      TRUE ~ NA_character_
    ),
    riparian_recommended = wetland_wettest %in% c("OBL", "FACW", "FAC")
  )

# =============================================================================
# STEP 4: CATERPILLAR HOST RICHNESS — NHM LONDON "HOSTS" DATABASE
# Confirmed column names from the actual downloaded CSV:
#   _id, Damage, Hostplant Family, Hostplant Genus, Hostplant Species,
#   Hostplant Subspecies/var, HOSTS ID, Insect Author, Insect Family,
#   Insect Genus, Insect Species, Insect Subspecies, Lab Rearing, Location
#
# Location values include: "USA", "Nearctic", "Canada", "Mexico",
#   "Puerto Rico", "North America", "British Isles", "Australia", etc.
# We keep USA + Canada + Mexico + Nearctic + Puerto Rico.
# "Nearctic" is the most important — it covers the entire NA bioregion
# and is how many literature records are coded regardless of exact country.
#
# Unique Lepidoptera species counted as distinct Insect Genus + Insect Species
# combinations per Hostplant Genus. Rows with blank Hostplant Genus are
# family-only records (mostly grasses/forbs) and are correctly excluded.
#
# NOTE on BAMONA (butterfliesandmoths.org): North America-focused but not
# usable here — bulk download requires a case-by-case data scientist account
# request and has redistribution restrictions. NHM HOSTS is CC0 and
# directly downloadable. Stick with HOSTS.
# Source: Robinson et al. (2023), CC0, Natural History Museum London
# Download: https://data.nhm.ac.uk/dataset/hosts  (HOSTS resource, CSV)
#
# The HOSTS CSV has columns including:
#   Hostplant genus, Hostplant species, Country, Lepidoptera family, etc.
# We filter to North American records and count distinct Lepidoptera species
# per hostplant genus — this reproduces the Tallamy methodology from a
# fully downloadable, citable, CC0 database.
# =============================================================================

if (file.exists(PATH_HOSTS_LEP)) {
  message("Loading NHM HOSTS database...")
  
  # Confirmed column names from the real HOSTS CSV (do not change):
  #   _id | Damage | Hostplant Family | Hostplant Genus | Hostplant Species |
  #   Hostplant Subspecies/var | HOSTS ID | Insect Author | Insect Family |
  #   Insect Genus | Insect Species | Insect Subspecies | Lab Rearing | Location
  
  hosts_raw <- read_csv(PATH_HOSTS_LEP, show_col_types = FALSE,
                        name_repair = "minimal")
  
  message(sprintf("HOSTS loaded: %d records", nrow(hosts_raw)))
  
  # Filter to North American records.
  # Location values in HOSTS: "USA", "Canada", "Mexico", "Nearctic",
  # "Puerto Rico", "North America", "Neotropical" (partial overlap), etc.
  # "Nearctic" is the biogeographic region covering all of North America
  # above Mexico and is used for many literature-sourced records — must include.
  na_locations <- c("usa", "canada", "mexico", "nearctic",
                    "north america", "puerto rico", "u.s.a",
                    "united states", "neotropical")  # neotropical catches Mexico/CA records
  
  hosts_na <- hosts_raw %>%
    filter(
      !is.na(Location),
      str_to_lower(Location) %in% na_locations |
        str_detect(str_to_lower(Location), "usa|nearctic|canada|mexico|puerto rico")
    )
  
  message(sprintf("HOSTS: %d North American records of %d total (%.0f%%)",
                  nrow(hosts_na), nrow(hosts_raw),
                  100 * nrow(hosts_na) / nrow(hosts_raw)))
  
  # Build unique Lepidoptera species key: Insect Genus + Insect Species.
  # Rows where Insect Species is blank are genus-only records — still count
  # as distinct Lep taxa at genus level if genus differs.
  hosts_na <- hosts_na %>%
    mutate(
      lep_key = paste(
        str_trim(replace_na(`Insect Genus`, "")),
        str_trim(replace_na(`Insect Species`, ""))
      ) %>% str_squish()
    )
  
  # Count distinct Lep species per hostplant genus.
  # Exclude rows where Hostplant Genus is blank (family-only plant records).
  lep_by_genus <- hosts_na %>%
    filter(
      !is.na(`Hostplant Genus`),
      str_trim(`Hostplant Genus`) != "",
      lep_key != ""
    ) %>%
    group_by(genus = str_trim(`Hostplant Genus`)) %>%
    summarise(
      lep_spp_supported = n_distinct(lep_key),
      lep_source = "NHM HOSTS DB (Robinson et al. 2023, CC0); North American Location records",
      .groups = "drop"
    ) %>%
    arrange(desc(lep_spp_supported))
  
  message(sprintf("Lep host richness computed for %d plant genera", nrow(lep_by_genus)))
  message("Top 10 genera by Lep richness:")
  print(head(lep_by_genus, 10))
  
} else {
  warning(paste(
    "HOSTS_lep.csv not found. Lepidoptera host data will be NA for all species.",
    "Download from: https://data.nhm.ac.uk/dataset/hosts",
    "Click the 'HOSTS' resource -> Download button -> save as HOSTS_lep.csv",
    "File is tab-separated despite the .csv extension.",
    sep = "\n"
  ))
  lep_by_genus <- tibble(genus = character(), lep_spp_supported = integer(),
                         lep_source = character())
}

# =============================================================================
# STEP 5: PEST / PATHOGEN FLAGS — USFS NIDRM MODEL DETAILS WORKBOOK
# Source: USDA Forest Service, public domain
# Download: https://www.fs.usda.gov/foresthealth/applied-sciences/
#             mapping-reporting/national-risk-maps
# Click: "2013-2027 NIDRM Model Details workbook" (Excel)
# Save as: NIDRM_pests.xlsx
#
# The workbook contains one row per pest:host pair with columns including
# pest name, host scientific name, damage type, and severity rating.
# We read it, filter to your tree species, and build the flag summary.
# =============================================================================

if (file.exists(PATH_NIDRM)) {
  message("Loading USFS NIDRM pest workbook...")
  
  # Sheet: "NIDRM Model Summary" — one row per pest-host model.
  # Confirmed columns (cleaned): model_id | model_name | host_code |
  #   dca_code | agent (pest common name) |
  #   agentcollapsed_for_3_groups (pest category) |
  #   host (host common name) | max_percent_mortality | region |
  #   model_certainty | bailey_eco_regions
  models <- read_excel(PATH_NIDRM, sheet = "NIDRM Model Summary",
                       col_types = "text") %>%
    clean_names() %>%
    rename(pest_name     = agent,
           pest_category = agen_tcollapsed_for_3_groups,
           max_mortality = max_percent_mortality) %>%
    mutate(max_mortality = as.numeric(max_mortality))
  
  message(sprintf("NIDRM Model Summary: %d models", nrow(models)))
  
  # Sheet: "Hosts" — confirmed columns: host_code | is_a_group_1_yes |
  #   host_common_name | abbreviation_used_in_model_name | genus | specific_epithet
  hosts_sheet <- read_excel(PATH_NIDRM, sheet = "Hosts",
                            col_types = "text") %>%
    clean_names() %>%
    rename(is_group  = is_a_group_1_yes,
           epithet   = specific_epithet) %>%
    mutate(
      scientific_name = if_else(
        is.na(is_group) & !is.na(genus) & !str_detect(coalesce(epithet,""), "^spp"),
        paste(str_trim(genus), str_trim(epithet)),
        NA_character_
      )
    )
  
  # Sheet: "Host Groups" — expands group host codes to individual species.
  # Confirmed columns: species_group_id | species_group_name |
  #   group_member_species_id | group_member_species_name
  host_groups <- read_excel(PATH_NIDRM, sheet = "Host Groups",
                            col_types = "text") %>%
    clean_names() %>%
    rename(group_code  = species_group_id,
           member_code = group_member_species_id)
  
  # Build flat host_code -> scientific_name lookup, including group expansion
  host_direct <- hosts_sheet %>%
    filter(is.na(is_group), !is.na(scientific_name)) %>%
    select(host_code, scientific_name)
  
  # Many group members (e.g. "Swamp White Oak", "Overcup Oak", "Chinkapin Oak")
  # never get their own row in the Hosts sheet -- NIDRM only gives individual
  # rows to species with their own dedicated pest models. The ID-based join
  # above (host_via_group) therefore misses most group members.
  #
  # Fix: build a common-name -> scientific-name lookup from our own species
  # list (species_df already has common_name from USDA PLANTS) and match
  # group_member_species_name against it directly. This recovers group
  # members that have no individual Hosts-sheet row.
  # common_name lives in plants_clean (built in STEP 2, joined to species_df
  # later in STEP 6) -- not in species_df itself at this point in the script.
  common_name_lookup <- plants_clean %>%
    filter(!is.na(common_name), !is.na(scientific_name)) %>%
    mutate(common_name_clean = str_to_lower(str_trim(common_name))) %>%
    distinct(common_name_clean, .keep_all = TRUE) %>%
    select(common_name_clean, scientific_name)
  
  host_via_group_by_id <- host_groups %>%
    left_join(host_direct, by = c("member_code" = "host_code")) %>%
    select(host_code = group_code, scientific_name) %>%
    filter(!is.na(scientific_name))
  
  host_via_group_by_name <- host_groups %>%
    mutate(common_name_clean = str_to_lower(str_trim(group_member_species_name))) %>%
    left_join(common_name_lookup, by = "common_name_clean") %>%
    select(host_code = group_code, scientific_name) %>%
    filter(!is.na(scientific_name))
  
  host_via_group <- bind_rows(host_via_group_by_id, host_via_group_by_name) %>%
    distinct(host_code, scientific_name)
  
  message(sprintf(
    "Group expansion: %d via host_code match, %d additional via common-name match",
    nrow(host_via_group_by_id),
    nrow(anti_join(host_via_group_by_name, host_via_group_by_id,
                   by = c("host_code","scientific_name")))
  ))
  
  host_lookup <- bind_rows(host_direct, host_via_group) %>%
    distinct(host_code, scientific_name)
  
  message(sprintf("Host lookup: %d host_code -> scientific_name pairs", nrow(host_lookup)))
  
  # Join models to scientific names, filter to our BA species
  pest_long <- models %>%
    left_join(host_lookup, by = "host_code", relationship = "many-to-many") %>%
    filter(!is.na(scientific_name))
  
  our_species <- species_df$scientific_name
  
  pest_raw <- pest_long %>%
    filter(scientific_name %in% our_species) %>%
    mutate(
      # Max % Mortality -> 1-4 severity scale
      severity = case_when(
        max_mortality >= 60 ~ 4L,
        max_mortality >= 30 ~ 3L,
        max_mortality >= 10 ~ 2L,
        TRUE                ~ 1L
      ),
      flag_user = severity >= 3,
      pest_type = pest_category,
      source    = "USFS NIDRM 2013-2027 Model Details Workbook"
    ) %>%
    select(scientific_name, pest_name, pest_type, max_mortality,
           severity, flag_user, region, source) %>%
    distinct()
  
  message(sprintf("NIDRM: %d pest-host records matched to your species", nrow(pest_raw)))
  if (nrow(pest_raw) > 0) {
    print(head(pest_raw %>% select(scientific_name, pest_name,
                                   max_mortality, severity, region), 10))
  }
  
} else {
  warning(paste(
    "NIDRM_pests.xlsx not found.",
    "Download from: https://www.fs.usda.gov/foresthealth/technology/docs/2012_NIDRM_Model_Info_12_04_14_external.xlsx",
    sep = "\n"
  ))
  pest_raw <- tibble(scientific_name = character(), pest_name = character(),
                     pest_type = character(), max_mortality = numeric(),
                     severity = integer(), flag_user = logical(),
                     region = character(), source = character())
}

# Collapse to one row per species for the main CSV
pest_summary <- pest_raw %>%
  group_by(scientific_name) %>%
  summarise(
    has_pest_flag     = any(flag_user, na.rm = TRUE),
    max_pest_severity = max(severity, na.rm = TRUE),
    # Flagged pests only in the summary list
    pest_flag_list    = paste(pest_name[flag_user], collapse = " | "),
    # Full detail: all pests with mortality % and region
    pest_flag_notes   = paste(
      paste0("[", pest_name, " — ", max_mortality, "% max mortality, sev.",
             severity, ", region:", coalesce(region, "?"), "]"),
      collapse = " | "),
    pest_source       = first(source),
    .groups = "drop"
  )

# =============================================================================
# STEP 6: JOIN EVERYTHING
# =============================================================================
join1 <- species_df %>%
  left_join(plants_clean,  by = "scientific_name")

join2 <- join1 %>%
  left_join(wetland_l48,   by = "scientific_name")

join3 <- join2 %>%
  left_join(lep_by_genus,  by = "genus")

join4 <- join3 %>%
  left_join(pest_summary,  by = "scientific_name") %>%
  mutate(
    has_pest_flag     = replace_na(has_pest_flag, FALSE),
    max_pest_severity = if_else(is.infinite(max_pest_severity), 0L,
                                as.integer(max_pest_severity)),
    pest_flag_list    = replace_na(pest_flag_list, "None in NIDRM"),
    pest_flag_notes   = replace_na(pest_flag_notes, "")
  )

# =============================================================================
# STEP 7: WILDLIFE VALUE INDEX (0-100)
# Identical logic to previous version — documented here for reference.
#
# A. MAST/FRUIT/SEED (0-25):
#    berry_nut_seed_product Yes=12  fruit_seed_abundance H=8/M=5/L=2
#    fruit_seed_persistence Yes=5   cap 25
#
# B. INSECT HOST RICHNESS (0-25):
#    lep_spp_supported / max_lep * 25  (max computed from data, not hardcoded)
#    NA -> imputed from mean of A, C, D
#
# C. BLOOM/POLLINATOR (0-25):
#    any bloom_period recorded = 15; spring = +5 bonus; cap 20
#    (pollen_or_nectar_source not in this DB version)
#
# D. VERTEBRATE FOOD/BROWSE (0-25):
#    palatable_browse_animal H=25/M=15/L=5; NA=12.5
#
# NATIVITY MULTIPLIER:
#    native L48 = 1.00 | non-native non-invasive = 0.85 | invasive = 0.60
# =============================================================================

# Max Lep value from the data itself (not hardcoded to Quercus 557)
max_lep <- if (nrow(lep_by_genus) > 0) max(lep_by_genus$lep_spp_supported, na.rm = TRUE) else 557

yn_pts <- function(f, y, n = 0) case_when(str_to_lower(trimws(f)) == "yes" ~ y,
                                          str_to_lower(trimws(f)) == "no"  ~ n,
                                          TRUE ~ NA_real_)
abund_pts  <- function(x) case_when(str_to_lower(trimws(x)) == "high" ~ 8,
                                    str_to_lower(trimws(x)) == "medium" ~ 5,
                                    str_to_lower(trimws(x)) == "low" ~ 2,
                                    TRUE ~ NA_real_)
browse_pts <- function(x) case_when(str_to_lower(trimws(x)) == "high" ~ 25,
                                    str_to_lower(trimws(x)) == "medium" ~ 15,
                                    str_to_lower(trimws(x)) == "low" ~ 5,
                                    TRUE ~ 12.5)
bloom_pts  <- function(x) {
  base  <- if_else(!is.na(x) & trimws(x) != "", 15, 0)
  bonus <- if_else(str_detect(tolower(x), "spring|early"), 5, 0, missing = 0)
  pmin(base + bonus, 20)
}

output_df <- join4 %>%
  mutate(
    cat_A = pmin(replace_na(yn_pts(berry_nut_seed_product, 12, 0), 0) +
                   replace_na(abund_pts(fruit_seed_abundance), 0) +
                   replace_na(yn_pts(fruit_seed_persistence, 5, 0), 0), 25),
    cat_B = pmin(lep_spp_supported / max_lep * 25, 25),
    cat_C = bloom_pts(bloom_period),
    cat_D = browse_pts(palatable_browse_animal),
    .b_fill   = if_else(is.na(cat_B), (cat_A + cat_C + cat_D) / 3, cat_B),
    raw_score = .b_fill + cat_A + cat_C + cat_D,
    is_l48_native  = replace_na(native_l48, FALSE),
    is_invasive_sp = str_detect(replace_na(pest_flag_list, ""), "INVASIVE|invasive"),
    nativity_mult  = case_when(is_invasive_sp ~ 0.60,
                               is_l48_native  ~ 1.00,
                               TRUE           ~ 0.85),
    wildlife_value_index = round(pmin(raw_score * nativity_mult, 100), 1)
  ) %>%
  select(
    ba_code, scientific_name, genus, common_name, is_unknown, name_verified,
    native_states, introduced_states, native_regions,
    # Wildlife value
    wildlife_value_index, cat_A, cat_B, cat_C, cat_D, nativity_mult, is_l48_native,
    lep_spp_supported, lep_source,
    berry_nut_seed_product, fruit_seed_abundance,
    fruit_seed_period_begin, fruit_seed_period_end, fruit_seed_persistence,
    palatable_browse_animal, bloom_period,
    # Riparian
    wetland_wettest, wetland_modal, riparian_category, riparian_recommended,
    all_of(avail_l48),
    # Ecology
    nitrogen_fixation, known_allelopath, anaerobic_tolerance,
    moisture_use, drought_tolerance, shade_tolerance, fire_tolerance,
    growth_rate, lifespan, resprout_ability, seedling_vigor,
    height_mature_ft, growth_form, growth_habits, active_growth_period, toxicity,
    # Soil
    adapted_to_coarse_textured_soils, adapted_to_fine_textured_soils,
    adapted_to_medium_textured_soils, soil_ph_min, soil_ph_max,
    salinity_tolerance, caco_3_tolerance,
    # Pest
    has_pest_flag, max_pest_severity, pest_flag_list, pest_flag_notes, pest_source
  ) %>%
  arrange(scientific_name)

# =============================================================================
# STEP 8: WRITE OUTPUT
# =============================================================================
write_csv(output_df, "BA_wildlife_traits.csv", na = "")
write_csv(pest_raw,  "BA_pest_pathogen_detail.csv", na = "")

message("\n=== DONE ===")
message(sprintf("Rows: %d | Columns: %d", nrow(output_df), ncol(output_df)))
message(sprintf("Tree filter applied — non-tree growth habits excluded"))
message(sprintf("Lepidoptera max value used for normalisation: %d", max_lep))
message(sprintf("Wildlife index range: %.1f – %.1f",
                min(output_df$wildlife_value_index, na.rm = TRUE),
                max(output_df$wildlife_value_index, na.rm = TRUE)))
message(sprintf("Riparian recommended: %d", sum(output_df$riparian_recommended, na.rm = TRUE)))
message(sprintf("Pest flagged: %d", sum(output_df$has_pest_flag)))
message(sprintf("Missing Lep data (cat_B NA): %d", sum(is.na(output_df$cat_B))))