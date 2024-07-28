BASE_URL:=http://hep.tsinghua.edu.cn/~berrysoft/bdeph2024/
WGET:=wget --quiet

geo.h5:
	$(WGET) $(BASE_URL)geo.h5

concat.h5:
	$(WGET) $(BASE_URL)concat.h5

draw.pdf: concat.h5
	./draw.py draw --concat $^ -o $@

.PHONY: score
score: concat.h5
	./draw.py validate --concat $^

seeds:=$(shell seq 16001 16020)

.PHONY: data
data: $(seeds:%=data/%.h5)

data/%.h5:
	@mkdir -p $(@D)
	$(WGET) -P $(@D) $(BASE_URL)$*.h5 
