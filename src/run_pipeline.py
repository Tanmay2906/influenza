"""Run stages 2-4 end to end (stage 1 needs network: python src/fetch_data.py)."""
import align_sequences, annotate_variants, call_variants, validate_sequences

for stage in (validate_sequences, align_sequences, call_variants, annotate_variants):
    print(f"== {stage.__name__}")
    stage.main()
