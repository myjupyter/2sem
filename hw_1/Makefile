all:
	./search.py -q 'На чикагской выставке 1893 года посетители катались на колесе обозрения, а на выставке 1933 года — на скайрайде (на илл.).' -l links.txt -n 10 -m 1 -o result1_1.json
	./search.py -q 'На чикагской выставке 1893 года посетители катались на колесе обозрения, а на выставке 1933 года — на скайрайде (на илл.).' -l links.txt -n 10 -m 0 -o result1_2.json
	./search.py -q 'На чикагской выставке 1893 года посетители катались на колесе обозрения, а на выставке 1933 года — на скайрайде (на илл.).' -l links.txt -n 10 -m 2 -o result1_3_lambda09.json
	./search.py -q 'На чикагской выставке 1893 года посетители катались на колесе обозрения, а на выставке 1933 года — на скайрайде (на илл.).' --lam 0.5 -l links.txt -n 10 -m 2 -o result1_3_lambda05.json
	./search.py -q 'В Средние века считалось, что если поэты не повторяют общеизвестные истины, то они лгут' -l links2.txt -n 10 -m 1 -o result2_1.json
	./search.py -q 'В Средние века считалось, что если поэты не повторяют общеизвестные истины, то они лгут' -l links2.txt -n 10 -m 0 -o result2_2.json
	./search.py -q 'В Средние века считалось, что если поэты не повторяют общеизвестные истины, то они лгут' --lam 0.5 -l links2.txt -n 10 -m 2 -o result2_3_lambda05.json
	./search.py -q 'В Средние века считалось, что если поэты не повторяют общеизвестные истины, то они лгут' -l links2.txt -n 10 -m 2 -o result2_3_lambda09.json
	./search.py -q 'Французские правые поломали стулья в кинотеатре, пытаясь сорвать показ реалистического фильма' --lam 0.5 -l links3.txt -n 10 -m 2 -o result3_3_lambda05.json
	./search.py -q 'Французские правые поломали стулья в кинотеатре, пытаясь сорвать показ реалистического фильма' -l links3.txt -n 10 -m 2 -o result3_3_lambda09.json
	./search.py -q 'Французские правые поломали стулья в кинотеатре, пытаясь сорвать показ реалистического фильма' -l links3.txt -n 10 -m 1 -o result3_1.json
	./search.py -q 'Французские правые поломали стулья в кинотеатре, пытаясь сорвать показ реалистического фильма' -l links3.txt -n 10 -m 0 -o result3_2.json	

ndcg:
	@echo "Запрос: Французские правые поломали стулья в кинотеатре, пытаясь сорвать показ реалистического фильма"
	@echo "tf.idf"
	@./ndcg.py -d templates/Сука_\(фильм,_1931\).json -r ./result3_1.json
	@echo "tf.idf"
	@./ndcg.py -d templates/Сука_\(фильм,_1931\).json -r ./result3_2.json
	@echo "prob. model (lambda=0.5)"
	@./ndcg.py -d templates/Сука_\(фильм,_1931\).json -r ./result3_3_lambda05.json
	@echo "prob. model (lambda=0.9)"
	@./ndcg.py -d templates/Сука_\(фильм,_1931\).json -r ./result3_3_lambda09.json
	@echo ""
	@echo "В Средние века считалось, что если поэты не повторяют общеизвестные истины, то они лгут"
	@echo "tf.idf"
	@./ndcg.py -d templates/Средние_века.json -d templates/Художественный_вымысел.json -r ./result2_1.json
	@echo "tf.idf"
	@./ndcg.py -d templates/Средние_века.json -d templates/Художественный_вымысел.json -r ./result2_2.json
	@echo "prob. model (lambda=0.5)"
	@./ndcg.py -d templates/Средние_века.json -d templates/Художественный_вымысел.json -r ./result2_3_lambda05.json
	@echo "prob. model (lambda=0.9)"
	@./ndcg.py -d templates/Средние_века.json -d templates/Художественный_вымысел.json -r ./result2_3_lambda09.json
	@echo ""
	@echo "На чикагской выставке 1893 года посетители катались на колесе обозрения, а на выставке 1933 года — на скайрайде (на илл.)."
	@echo "tf.idf"
	@./ndcg.py -d templates/Всемирная_выставка_\(1893\).json -d templates/Всемирная_выставка_\(1933\).json -d templates/Колесо_обозрения.json -d templates/Скайрайд.json -r ./result1_1.json
	@echo "tf.idf"
	@./ndcg.py -d templates/Всемирная_выставка_\(1893\).json -d templates/Всемирная_выставка_\(1933\).json -d templates/Колесо_обозрения.json -d templates/Скайрайд.json -r ./result1_2.json
	@echo "prob. model (lambda=0.5)"
	@./ndcg.py -d templates/Всемирная_выставка_\(1893\).json -d templates/Всемирная_выставка_\(1933\).json -d templates/Колесо_обозрения.json -d templates/Скайрайд.json -r ./result1_3_lambda05.json
	@echo "prob. model (lambda=0.9)"
	@./ndcg.py -d templates/Всемирная_выставка_\(1893\).json -d templates/Всемирная_выставка_\(1933\).json -d templates/Колесо_обозрения.json -d templates/Скайрайд.json -r ./result1_3_lambda09.json
