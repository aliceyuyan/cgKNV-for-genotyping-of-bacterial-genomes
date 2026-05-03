/home/yichen/alice/software/matlab2020b/bin/matlab  --nosplash   --nodesktop

clc
clear
mkdir kmer

core_gene_position = 'Core_gene_positions.txt';
txtc1 ='';
txtc1 = textread(core_gene_position,'%s','delimiter','\n','whitespace','','bufsize', 4095000);
Starts=[];
Ends=[];
for x2=2:length(txtc1)
    content='';
    content=txtc1{x2};
    sites=[];
    sites=strfind(content,9);
    Starts(x2-1)=str2num(content(sites(1)+1:sites(2)-1));
    Ends(x2-1)=str2num(content(sites(2)+1:sites(3)-1));
end
clearvars -except Starts Ends

txtc2 ='';
txtc2 = textread('../Dataset2_list_583.txt','%s','delimiter','\n','whitespace','','bufsize', 4095000);
name_list='';
for x=1:length(txtc2)%%3083
    content='';
    content=txtc2{x};
    name_list{x}=content;
end

for i=1:length(name_list)
    i
    data_name='';
    data_name=name_list{i};
    consensus_file='';
    consensus_file=strcat(data_name,'.consensus.subs.fa');
    txtc = textread(consensus_file,'%s','delimiter','\n','whitespace','','bufsize', 4095000);
    mt = find(cellfun('isempty',txtc));
    txtc(mt) = [];
    lengthn = size(txtc, 1);
    taxanum = 0;
    for x = 1 : lengthn
        if txtc{x}(1) == '>'
            taxanum = taxanum + 1;
        end
    end
    seq='';
    seq{taxanum} ='';
    signk = 0;
    signtaxa = 1;
    for x1 = 2 : lengthn
        if txtc{x1}(1) ~= '>'
            lengthj = length(txtc{x1});
            seq{signtaxa}(signk+1 : signk+lengthj) = txtc{x1};
            signk = signk + lengthj;
        else
            signtaxa = signtaxa + 1;
            signk = 0;
        end
    end
    Consensus='';
    Consensus=seq{1};
    length(Consensus)%5227419

    S=[];
    for j=1:length(Starts)
        seqs='';
        seqs=upper(Consensus(Starts(j):Ends(j)));
        kmer=ceil(log(length(seqs))/log(4));
        nucleotide='ACGT';
        base=nucleotide';
        A='';
        for k=1:(kmer-1)
            l=1;
            for m=1:4^k
                for n=1:4
                    A(l,:)=strcat(base(m,:),nucleotide(n));
                    l=l+1;
                end
            end
            base=A;
            A='';
        end
        B=base;
        Nk=4^kmer;

        len=length(seqs)-kmer+1;
        F=[];
        for k=1:Nk
            C='';
            C=B(k,:);
            D=[];
            D=strfind(seqs,C);

            E=[];
            if length(D)~=0
                E=[length(D) mean(D) sum((D-mean(D)).^2)/(length(D)*len)];
            else
                E=[0 0 0];
            end
            F=[F E];
        end

        S=[S F];
    end

    save(strcat('./kmer/',data_name,'_kmer.mat'),'S');

end

cd kmer
measures={'hamming'};

for i1=1:length(measures)

    measure='';
    measure=measures{i1};

    V=[];
    for i=1:length(name_list)-1
        for j=i+1:length(name_list)
            data_name1='';
            data_name1=strcat(name_list{i},'_kmer.mat');
            load(data_name1);
            m=[];
            m=S;
            clearvars S

            data_name2='';
            data_name2=strcat(name_list{j},'_kmer.mat');
            load(data_name2);
            n=[];
            n=S;
            clearvars S

            V(i,j)=pdist([m;n],measure);
            clearvars m n
        end
    end

    taxanum=length(name_list);
    M=zeros(taxanum);
    for i=1:taxanum-1
        for j=i+1:taxanum
            M(i,j)=V(i,j);
            M(j,i)=M(i,j);
        end
    end

    save(strcat('anthracis_cgKNV_trimmomatic_',measure,'.mat'),'M');

    outputfilename=strcat('anthracis_cgKNV_trimmomatic_',measure,'.meg');
    num_of_seq=taxanum;
    outfile=fopen(outputfilename, 'w');
    fprintf(outfile, '#mega\n');
    fprintf(outfile, '!Title: TEST;\n');
    fprintf(outfile, '!Format DataType=Distance DataFormat=LowerLeft NTaxa=%d;\n', num_of_seq);
    fprintf(outfile, '\n');
    for k = 1 : num_of_seq
        fprintf(outfile, '[%d] #%s\n', k, name_list{k});
    end
    fprintf(outfile, '\n');
    for j = 2 : num_of_seq
        fprintf(outfile, '[%d]   ', j);
        for k = 1 : (j-1)
            fprintf(outfile, ' %8f', M(j, k));
        end
        fprintf(outfile, '\n');
    end
    fprintf(outfile, '\n');
    fclose(outfile);

end

