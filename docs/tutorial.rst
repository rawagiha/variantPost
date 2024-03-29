.. _Tutorial:

Tutorial
=========


Querying VCF file
-----------------

|

Using del2 (bottom) as example, query the COSMIC VCF database::
    
    import pysam
    from variantpost import Variant
     
    reference = pysam.FastaFile("/path/to/GRCh38.fa")
    cosmic = pysam.VariantFile("/path/to/cosmic.v89.vcf(.gz)")

    # del2 
    v = Variant("17", 31224665, "CC", "C", reference)
    
Normalization query (default) returns VCF entries that are identical after normalization::
    
    norm_hits = v.query_vcf(cosmic) # list of 2 hit VCF entries (del1 and del2)
    
    for hit in norm_hits:
        print(hit["INFO"]["CNT"]) 
        
        #COSMIC count for del1
        #COSMIC count for del2  

Locus query returns VCF entries located at the normalized genomic locus::

    locus_hits = v.query_vcf(cosmic, matchby="locus") # list of 5 hit VCF entries (all indels)

    for hit in locus_hits:
        print(hit["INFO"]["CNT"]) 
        
        #COSMIC count for del1
        ...
        #COSMIC count for ins3

Exact query only returns a VCF entry matching without normalization:: 
        
    exact_hit = v.query(cosmic, matchby="exact") # list of a hit VCF entry (del2)
    
    print(exact_hit[0]["INFO"]["CNT"]) 
    
    #COSMIC count for del2
    

Phasing to complex variant
---------------------------
