import igv
b = igv.Browser({"genome": "hg38"})
b.show()

b.load_track(
    {
        # "name": "Cerebral cortex pl",
        # "url": "data/alignmentmerge_CHTN22a_and_CHTN22b_pl.bw",
        "name": "Astrocytes pl",
        "url": "data/alignmentmerge_SC3a_and_SC3b_pl.bw",
        # "name": "Neurons pl",
        # "url": "data/alignmentmerge_SC2a_and_SC2b_pl.bw",
        "type": "wig",
        "indexed": False,
        "color": "red"
    })
b.load_track(
    {
        # "name": "Cerebral cortex mn",
        # "url": "data/alignmentmerge_CHTN22a_and_CHTN22b_mn.bw",
        "name": "Astrocytes mn",
        "url": "data/alignmentmerge_SC3a_and_SC3b_mn.bw",
        # "name": "Neurons mn",
        # "url": "data/alignmentmerge_SC2a_and_SC2b_mn.bw",
        "type": "wig",
        "indexed": False,
        "color": "blue"
    })
b.load_track(
    {
        "name": "Enhancer elements",
        "url": "data/linked_enhancers.bed",
        "type": "annotation",
        "indexed": False,
        "color": "#2e57ac"
    })
b.load_track(
    {
        "name": "DNV_case",
        "url": "data/enhancer_DNV_case.bed",
        "type": "annotation",
        "indexed": False,
        "color": "red"
    })
b.load_track(
    {
        "name": "DNV_control",
        "url": "data/enhancer_DNV_control.bed",
        "type": "annotation",
        "indexed": False,
        "color": "green"
    })
b.load_track(
    {
        "name": "Transmitted_case",
        "url": "data/enhancer_Transmitted_case.bed",
        "type": "annotation",
        "indexed": False,
        "color": "red"
    })
b.load_track(
    {
        "name": "Untransmitted_case",
        "url": "data/enhancer_Untransmitted_case.bed",
        "type": "annotation",
        "indexed": False,
        "color": "green"
    })
# chr1:58773306-58773839 # "chr1:58772806-58774939" (cerebral cortex)
# chr1:58779786-58780532 # "chr1:58779386-58780532" (astrocytes)
# chr1:58780876-58781490 # "chr1:58780336-58781490" (astrocytes)
# chr1:58781631-58783539 # "chr1:58778631-58783939" (cerebral cortex)
# chr1:58784641-58786158 # "chr1:58784450-58786158" (cerebral cortex)
b.search("chr1:58780336-58781490")