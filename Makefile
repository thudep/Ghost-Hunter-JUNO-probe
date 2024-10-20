BASE_URL:=https://ghfile.thudep.com:4433
WGET:=wget

geo.h5:
	$(WGET) $(BASE_URL)/geo/geo.h5

concat.h5:
	$(WGET) $(BASE_URL)/test/concat.h5

draw.pdf: concat.h5 histogram.h5
	python3 draw.py draw --concat $< -o $@

.PHONY: score
score: concat.h5 histogram.h5
	python3 draw.py validate --concat $<

histogram.h5: geo.h5 data
	python3 histogram.py -g $< --data $(word 2,$^) -o $@ -b 10 -t 10

seeds:=$(shell seq 16001 16001)

.PHONY: data
data: $(seeds:%=data/%.h5)

clean:
	rm -rf data/
	rm -rf *.pdf
	rm -rf *.h5
	rm -rf __pycache__

data/%.h5:
	@mkdir -p $(@D)
	$(WGET) -P $(@D) $(BASE_URL)/data/$*.h5
