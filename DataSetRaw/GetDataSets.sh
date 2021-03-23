mkdir web-Stanford
mkdir dblp
mkdir pokec
mkdir liveJournal
mkdir orkut
mkdir twitter
cd web-Stanford     &&  wget https://snap.stanford.edu/data/web-Stanford.txt.gz                             && cd ..
cd dblp             &&  wget https://snap.stanford.edu/data/bigdata/communities/com-dblp.ungraph.txt.gz     && cd ..
cd pokec            &&  wget https://snap.stanford.edu/data/soc-pokec-relationships.txt.gz                  && cd ..
cd liveJournal      &&  wget https://snap.stanford.edu/data/soc-LiveJournal1.txt.gz                         && cd ..
cd orkut            &&  wget https://snap.stanford.edu/data/bigdata/communities/com-orkut.ungraph.txt.gz    && cd ..
cd twitter          &&  wget https://snap.stanford.edu/data/twitter-2010.txt.gz                             && cd ..


