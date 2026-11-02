import importlib.util, json, os, tempfile, unittest
from pathlib import Path
from unittest.mock import patch

s = importlib.util.spec_from_file_location('worker', Path(__file__).with_name('mapping_worker.py'))
w = importlib.util.module_from_spec(s); s.loader.exec_module(w)

class GenomeOnlyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.root = Path(self.temp.name)
        self.work = self.root/'work'; self.work.mkdir()
        self.genome = self.root/'genome.fa'; self.genome.write_text('>chrX\nACGTACGT\n')
        self.c = dict(assembly='GCA_004329205.1_UCBerk_Dpse_1.0',
            analysis_mode='genome_mapping_only', annotation_status='No matching annotation',
            genome=str(self.genome), genome_sha256=w.digest(self.genome),
            reference_prefix='Dpse', reference_dir=str(self.root/'reference'),
            results_dir=str(self.root/'results'))
        self.tools = {'STAR':{'version':'STAR-test'}}
    def tearDown(self): self.temp.cleanup()

    def test_reference_build_requires_no_annotation_or_rsem(self):
        commands = []
        def fake_run(args, work, log, stdout=None):
            commands.append(list(map(str,args)))
            self.assertEqual(args[0], 'STAR')
            for name in ('Genome','SA','SAindex','genomeParameters.txt'):
                (work/name).write_text('test-index')
        with patch.dict(os.environ, SLURM_JOB_ID='toy-reference'), patch.object(w,'run',fake_run):
            w.build_reference(self.c,self.work,8,self.tools)
        self.assertNotIn('--sjdbGTFfile',commands[0])
        receipt=w.reference_check(self.c,with_hashes=True)
        self.assertEqual(receipt['assembly'],self.c['assembly'])
        self.assertEqual(set(receipt['source_sha256']),{'genome'})
        self.c['genome_sha256']='wrong'
        with self.assertRaises(ValueError): w.reference_check(self.c)

    def test_mapping_outputs_bam_without_claiming_gene_quantification(self):
        import gzip
        paths=[]
        for mate in (1,2):
            p=self.root/('r%d.gz'%mate)
            with gzip.open(p,'wb') as f:f.write(b'@one\nACGT\n+\nIIII\n')
            paths.append(str(p))
        sample=dict(sample_id='Dpse.test', run_accession='SRRtest',layout='PE',
                    strandedness='auto',read1=paths[0],read2=paths[1],gzip_depths=[1,1],
                    input_metadata=[w.metadata(p) for p in paths])
        commands=[]
        def fake_run(args,work,log,stdout=None):
            args=list(map(str,args)); commands.append(args)
            self.assertFalse(args[0].startswith('rsem'))
            if args[0]=='STAR':
                (work/'star.Log.final.out').write_text('Number of input reads | 1\n')
                (work/'star.Aligned.out.bam').write_bytes(b'toy-bam')
                self.assertNotIn('--quantMode',args)
                self.assertEqual(args[args.index('--twopassMode')+1],'Basic')
            elif args[1]=='sort':Path(args[args.index('-o')+1]).write_bytes(b'toy-sorted-bam')
            elif args[1]=='index':Path(args[-1]+'.bai').write_bytes(b'toy-index')
            elif args[1]=='quickcheck':
                for p in args[3:]:self.assertTrue(Path(p).exists(),p)
            if stdout:Path(stdout).write_text('2\n' if args[1]=='view' else 'toy-QC\n')
        with patch.dict(os.environ,SLURM_JOB_ID='toy-sample'), \
             patch.object(w,'reference_check',return_value={'tools':self.tools}), \
             patch.object(w,'run',fake_run), \
             patch.object(w,'infer_strand',side_effect=AssertionError('No annotation inference permitted')):
            self.assertEqual(w.map_sample(self.c,sample,self.work,8,self.tools),0)
        result=self.root/'results'/'Dpse.test'
        receipt=json.loads((result/'receipt.json').read_text())
        self.assertEqual(receipt['state'],'GENOME_MAPPED_NO_GENE_QUANTIFICATION')
        self.assertTrue((result/'Dpse.test.genome.bam').exists())
        self.assertTrue((result/'Dpse.test.genome.bam.bai').exists())
        self.assertFalse(list(result.glob('*.genes.results')))
        self.assertEqual(receipt['fastq_integrity_check']['checked_records_per_file'],1)

if __name__=='__main__':unittest.main(verbosity=2)
