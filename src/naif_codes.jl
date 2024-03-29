# From https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
# TODO do we need the rest?
export naif_body_codes, naif_reference_codes
naif_body_codes = Dict(
                       # Barycenters
                       "SOLAR_SYSTEM_BARYCENTER" => 0, "SSB" => 0,
                       "SOLAR SYSTEM BARYCENTER" => 0, "MERCURY_BARYCENTER" => 1,
                       "MERCURY BARYCENTER" => 1, "VENUS_BARYCENTER" => 2,
                       "VENUS BARYCENTER" => 2, "EARTH_BARYCENTER" => 3, "EMB" => 3,
                       "EARTH MOON BARYCENTER" => 3, "EARTH-MOON BARYCENTER" => 3,
                       "EARTH BARYCENTER" => 3, "MARS_BARYCENTER" => 4,
                       "MARS BARYCENTER" => 4, "JUPITER_BARYCENTER" => 5,
                       "JUPITER BARYCENTER" => 5, "SATURN_BARYCENTER" => 6,
                       "SATURN BARYCENTER" => 6, "URANUS_BARYCENTER" => 7,
                       "URANUS BARYCENTER" => 7, "NEPTUNE_BARYCENTER" => 8,
                       "NEPTUNE BARYCENTER" => 8, "PLUTO_BARYCENTER" => 9,
                       "PLUTO BARYCENTER" => 9, "SUN" => 10,
                       # Planets and Satellites
                       "MERCURY" => 199, "VENUS" => 299, "EARTH" => 399, "MOON" => 301,
                       "MARS" => 499, "PHOBOS" => 401, "DEIMOS" => 402, "JUPITER" => 599,
                       "IO" => 501, "EUROPA" => 502, "GANYMEDE" => 503, "CALLISTO" => 504,
                       "AMALTHEA" => 505, "HIMALIA" => 506, "ELARA" => 507,
                       "PASIPHAE" => 508, "SINOPE" => 509, "LYSITHEA" => 510,
                       "CARME" => 511, "ANANKE" => 512, "LEDA" => 513, "THEBE" => 514,
                       "ADRASTEA" => 515, "METIS" => 516, "CALLIRRHOE" => 517,
                       "THEMISTO" => 518, "MAGACLITE" => 519, "TAYGETE" => 520,
                       "CHALDENE" => 521, "HARPALYKE" => 522, "KALYKE" => 523,
                       "IOCASTE" => 524, "ERINOME" => 525, "ISONOE" => 526,
                       "PRAXIDIKE" => 527, "AUTONOE" => 528, "THYONE" => 529,
                       "HERMIPPE" => 530, "AITNE" => 531, "EURYDOME" => 532,
                       "EUANTHE" => 533, "EUPORIE" => 534, "ORTHOSIE" => 535,
                       "SPONDE" => 536, "KALE" => 537, "PASITHEE" => 538, "HEGEMONE" => 539,
                       "MNEME" => 540, "AOEDE" => 541, "THELXINOE" => 542, "ARCHE" => 543,
                       "KALLICHORE" => 544, "HELIKE" => 545, "CARPO" => 546,
                       "EUKELADE" => 547, "CYLLENE" => 548, "KORE" => 549, "HERSE" => 550,
                       "DIA" => 553, "SATURN" => 699, "MIMAS" => 601, "ENCELADUS" => 602,
                       "TETHYS" => 603, "DIONE" => 604, "RHEA" => 605, "TITAN" => 606,
                       "HYPERION" => 607, "IAPETUS" => 608, "PHOEBE" => 609, "JANUS" => 610,
                       "EPIMETHEUS" => 611, "HELENE" => 612, "TELESTO" => 613,
                       "CALYPSO" => 614, "ATLAS" => 615, "PROMETHEUS" => 616,
                       "PANDORA" => 617, "PAN" => 618, "YMIR" => 619, "PAALIAQ" => 620,
                       "TARVOS" => 621, "IJIRAQ" => 622, "SUTTUNGR" => 623, "KIVIUQ" => 624,
                       "MUNDILFARI" => 625, "ALBIORIX" => 626, "SKATHI" => 627,
                       "ERRIAPUS" => 628, "SIARNAQ" => 629, "THRYMR" => 630, "NARVI" => 631,
                       "METHONE" => 632, "PALLENE" => 633, "POLYDEUCES" => 634,
                       "DAPHNIS" => 635, "AEGIR" => 636, "BEBHIONN" => 637,
                       "BERGELMIR" => 638, "BESTLA" => 639, "FARBAUTI" => 640,
                       "FENRIR" => 641, "FORNJOT" => 642, "HATI" => 643, "HYRROKKIN" => 644,
                       "KARI" => 645, "LOGE" => 646, "SKOLL" => 647, "SURTUR" => 648,
                       "ANTHE" => 649, "JARNSAXA" => 650, "GREIP" => 651, "TARQEQ" => 652,
                       "AEGAEON" => 653, "URANUS" => 799, "ARIEL" => 701, "UMBRIEL" => 702,
                       "TITANIA" => 703, "OBERON" => 704, "MIRANDA" => 705,
                       "CORDELIA" => 706, "OPHELIA" => 707, "BIANCA" => 708,
                       "CRESSIDA" => 709, "DESDEMONA" => 710, "JULIET" => 711,
                       "PORTIA" => 712, "ROSALIND" => 713, "BELINDA" => 714, "PUCK" => 715,
                       "CALIBAN" => 716, "SYCORAX" => 717, "PROSPERO" => 718,
                       "SETEBOS" => 719, "STEPHANO" => 720, "TRINCULO" => 721,
                       "FRANCISCO" => 722, "MARGARET" => 723, "FERDINAND" => 724,
                       "PERDITA" => 725, "MAB" => 726, "CUPID" => 727, "NEPTUNE" => 899,
                       "TRITON" => 801, "NEREID" => 802, "NAIAD" => 803, "THALASSA" => 804,
                       "DESPINA" => 805, "GALATEA" => 806, "LARISSA" => 807,
                       "PROTEUS" => 808, "HALIMEDE" => 809, "PSAMATHE" => 810, "SAO" => 811,
                       "LAOMEDEIA" => 812, "NESO" => 813, "PLUTO" => 999, "CHARON" => 901,
                       "NIX" => 902, "HYDRA" => 903, "KERBEROS" => 904, "STYX" => 905)

# https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/frames.html
naif_reference_codes = Dict("J2000" => 1, "B1950" => 2, "FK4" => 3, "DE-118" => 4,
                            "DE-96" => 5, "DE-102" => 6, "DE-108" => 7, "DE-111" => 8,
                            "DE-114" => 9, "DE-122" => 10, "DE-125" => 11, "DE-130" => 12,
                            "GALACTIC" => 13, "DE-200" => 14, "DE-202" => 15,
                            "MARSIAU" => 16, "ECLIPJ2000" => 17, "ECLIPB1950" => 18,
                            "DE-140" => 19, "DE-142" => 20, "DE-143" => 21)
