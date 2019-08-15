python -c "import logmuse
print(logmuse.LEVEL_BY_VERBOSITY)"





parser = logmuse.add_logmuse_args(parser)
args = parser.parse_args()

lmargs = logmuse.retrieve_logmuse_args(args)


global _LOGGER 
_LOGGER = logmuse.init_logger(lmargs)
logmuse.init_logger(name="pararead", lmargs)



non-CLI packages like pypiper and pararead should not need to know about logmuse,
but when a CLI tool (like a pipeline, or a count_reads tool) imports one of these with a logger, and uses logmuse, it should "just work".




file = "../microtest/data/bs_aln_k1k3.bam"