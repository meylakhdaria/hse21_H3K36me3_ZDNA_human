# hse21_H3K36me3_ZDNA_human

*Выполнила Мейлах Дарья*

## Анализ пиков гистоновой метки

### Выбор экспериментов

Для анализа я выбрала гистоновую метку **H3K36me3** с типом клеток **H9** организма **human (hg38 с приведением к hg19)** и вторичной структуры **ZDNA_DeepZ**.

Я выбрала два .bed файла ChIP-seq экспериментов из ENCODE: *ENCFF398TIJ* и *ENCFF905GSB*.

Скачала их с помощью команд.

*wget https://www.encodeproject.org/files/ENCFF398TIJ/@@download/ENCFF398TIJ.bed.gz*

*wget https://www.encodeproject.org/files/ENCFF905GSB/@@download/ENCFF905GSB.bed.gz*

### Необходимая обработка

И обрезала необходимое количество столбцов.

*zcat ENCFF398TIJ.bed.gz | cut -f1-5 > H3K36me3.ENCFF398TIJ.hg38.bed*

*zcat ENCFF905GSB.bed.gz | cut -f1-5 > H3K36me3.ENCFF905GSB.hg38.bed*

Далее провела конвертацию координат ChIP-seq пиков к версии генома hg19.

*liftOver   H3K36me3.ENCFF398TIJ.hg38.bed   hg38ToHg19.over.chain.gz   H3K36me3.ENCFF398TIJ.hg19.bed   H3K36me3.ENCFF398TIJ.unmapped.bed*

*liftOver   H3K36me3.ENCFF905GSB.hg38.bed   hg38ToHg19.over.chain.gz   H3K36me3.ENCFF905GSB.hg19.bed   H3K36me3.ENCFF905GSB.unmapped.bed*

### Построение гистограмм

Строю гистограмму длин участков для каждого эксперимента до и после конвертации с помощью [скрипта](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/src/len_hist.R).

**ENCFF398TIJ**:

[hg19](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/results/len_hist.H3K36me3.ENCFF398TIJ.hg19.pdf) - 35162 пика.

[hg38](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/results/len_hist.H3K36me3.ENCFF398TIJ.hg38.pdf) - 35170 пика.


**ENCFF905GSB**:

[hg19](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/results/len_hist.H3K36me3.ENCFF905GSB.hg19.pdf) - 15417 пика.

[hg38](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/results/len_hist.H3K36me3.ENCFF905GSB.hg38.pdf) - 15424 пика.

### Фильтрация

Среди ChIP-seq пиков для hg19  выкидываю пики >= 50000 с помощью [скрипта](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/src/filter.R) и строю для полученных данных гистограмму длин участков.

[ENCFF398TIJ hg19 filtered](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/results/filter_peaks.H3K36me3.ENCFF398TIJ.hg19.filtered.hist.pdf) - 35099 пика.

[ENCFF905GSB hg19 filtered](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/results/len_hist.H3K36me3.ENCFF905GSB.hg19.filtered.pdf) - 15417 пика.

Смотрю, где располагаются пики гистоновой метки относительно аннотированных генов с помощью [скрипта](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/src/pie.R).

**ENCFF398TIJ filtered**
![ENCFF398TIJ hg19 filtered](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/results/chip_seeker.H3K36me3.ENCFF398TIJ.hg19.filtered.plotAnnoPie.png)

**ENCFF905GSB filtered**
![ENCFF905GSB hg19 filtered](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/results/chip_seeker.H3K36me3.ENCFF905GSB.hg19.filtered.plotAnnoPie.png)

### Объединение

Объединяем два набора отфильтрованных ChIP-seq пиков с помощью утилиты bedtools merge.

*cat  *.filtered.bed  |   sort -k1,1 -k2,2n   |   bedtools merge   >  H3K36me3.merge.hg19.bed*

Визуализирую исходные два набора ChIP-seq пиков, а также их объединение в геномном браузере, и проверяю корректность работы bedtools merge.

Результат корректен. 
![gb](http://genome.ucsc.edu/trash/hgt/hgt_genome_50c4a_93970.png)

*track visibility=dense name="ENCFF398TIJ"  description="H3K36me3.ENCFF398TIJ.hg19.filtered.bed"
https://raw.githubusercontent.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/main/data/H3K36me3.ENCFF398TIJ.hg19.filtered.bed*

*track visibility=dense name="ENCFF905GSB"  description="H3K36me3.ENCFF905GSB.hg19.filtered.bed"
https://raw.githubusercontent.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/main/data/H3K36me3.ENCFF905GSB.hg19.filtered.bed*

*track visibility=dense name="ChIP_merge"  color=50,50,200   description="H3K36me3.merge.hg19.bed"
https://raw.githubusercontent.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/main/data/H3K36me3.merge.hg19.bed*

## Анализ участков вторичной структуры ДНК

Скачиваю файл со вторичной структурой ДНК ZDNA_DeepZ.

*wget https://raw.githubusercontent.com/Nazar1997/DeepZ/master/annotation/DeepZ.bed*

Строю распределение длин участков вторичной структуры ДНК с помощью [скрипта](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/src/len_hist.R).

[ZDNA_DeepZ](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/results/len_hist.DeepZ.pdf) - 19394 пиков.

Смотрю, где располагаются участки структуры ДНК относительно аннотированных генов с помощью [скрипта](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/src/pie.R).

![zdna](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/results/chip_seeker.DeepZ.plotAnnoPie.png)

## Анализ пересечений гистоновой метки и структуры ДНК

### Пересечения

Нахожу пересечения гистоновой меткой и структурой ДНК ZDNA_DeepZ

bedtools intersect  -a DeepZ.bed   -b  H3K36me3.merge.hg19.bed  >  H3K36me3.intersect_with_DeepZ.bed

Строю распределение длин с помощью [скрипта](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/src/len_hist.R).

[H3K36me3 intersect](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/results/len_hist.H3K6me36.intersect_with_DeepZ.pdf) - 147 пиков.

Визуализирую все в геномном браузере. Для этого добавляю для отображения следующие команды.

*track visibility=dense name="DeepZ"  color=0,200,0  description="DeepZ"
https://raw.githubusercontent.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/main/data/DeepZ.bed*

*track visibility=dense name="intersect_with_DeepZ"  color=200,0,0  description="H3K6me36.intersect_with_DeepZ.bed"
https://raw.githubusercontent.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/main/data/H3K6me36.intersect_with_DeepZ.bed*

Теперь мы можем отследить пересечения структуры и гистоновой метки, и так же их расположение относительно аннтированных генов. Вся сессия сохранена [здесь](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A23884564-23885763&hgsid=1124025283_XqXHaWSQBPUXiB2XY8RMb9Lf2bi9).

Приведу некоторые примеры пересечений.

![1](http://genome.ucsc.edu/trash/hgt/hgt_genome_33ab3_acc30.png)

![2](http://genome.ucsc.edu/trash/hgt/hgt_genome_3169b_abdc0.png)

### Ассоциация пересечений

Ассоциирую полученные пересечения с ближайшими генами с помощью  [скрипта](https://github.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/blob/main/src/anno.R).

Получаю [файл ассоциаций пиков с генами](https://raw.githubusercontent.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/main/data/H3K6me36.intersect_with_DeepZ.genes.txt), а также [список уникальных генов](https://raw.githubusercontent.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/main/data/H3K6me36.intersect_with_DeepZ.genes_uniq.txt). 

Всего было проассоциирован **41** пик с генами при общем количества уникальных генов равном **32**.

### GO-анализ

Провожу GO-анализ на http://pantherdb.org/
Полный список результатов представлен [здесь](https://raw.githubusercontent.com/meylakhdaria/hse21_H3K36me3_ZDNA_human/main/data/pantherdb_GO_analysis.txt).

Наиболее значимыми оказались категории  следующими значениями FDR:

*cotranslational protein targeting to membrane	FDR=7.69E-08*

*establishment of protein localization to endoplasmic reticulum	FDR=7.77E-08*

*cytoplasmic translation	FDR=7.78E-08*
