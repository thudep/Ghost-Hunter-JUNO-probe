concat.h5:
	# TODO
	# wget xxx

draw.pdf: concat.h5
	./draw.py draw --concat $^ -o $@

.PHONY: score
score: concat.h5
	@./draw.py validate --concat $^
