use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $reg = "Bio::EnsEMBL::Registry";

$reg->load_registry_from_db(
    -host   => 'ensembldb.ensembl.org',
    -user   => 'anonymous',
    -dbname => 'homo_sapiens_core_75_37'
);

open(OUTFILE, ">human_canonical_transcripts_v75_37");

my $gene_adaptor = $reg->get_adaptor( 'homo_sapiens', 'core', 'gene' );
my @genes = @{$gene_adaptor->fetch_all};
my $count = 0;

while(my $gene = shift @genes) {
    print
        ( OUTFILE
        $gene->external_name(), ":", $gene->stable_id, ":",$gene->canonical_transcript->stable_id, ".",
        $gene->canonical_transcript->version, "\n" );
}

close (OUTFILE);