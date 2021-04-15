      data LMHMbs(5) /11.7000/

      data (RGBb_coef(5,i),i=1,nRGBbc) /
     &     +2.406370D-01,
     &     +1.089220D+00,
     &     +1.953180D+00,
     &     -1.032180D+00,
     &     +0.000000D+00/

!     0.02000     LMR1
      data (mrange(5,1,i),i=1,3) /0,1,1/
      data (rrange(5,1,i),i=1,3) /0,3,1/
      data (alphas(5,1,i),i=1,8) /
     &     +1.50305864307172658556D+01,
     &     +4.96110465703095582235D-01,
     &     -9.16608282566077403608D-01,
     &     +2.48033316991016244968D-01,
     &     +1.77238406649739155263D+00,
     &     -6.52669943496246296455D-01,
     &     +6.23616025599144307989D-01,
     &     -1.77703681033998805994D-01/
      
!     0.02000     LMR2
      data (mrange(5,2,i),i=1,3) /0,5,1/
      data (rrange(5,2,i),i=1,3) /0,5,1/
      data (alphas(5,2,i),i=1,36) /
     &     +1.56617552673034907684D+01,
     &     -3.38271409511212528543D+00,
     &     +3.99712000854973803499D+00,
     &     -3.02506056136319223526D+00,
     &     +1.12651646333780397491D+00,
     &     -1.62914377442334834534D-01,
     &     +7.36910246144982838956D+00,
     &     -2.30880577686855161801D+01,
     &     +3.32294193400783512971D+01,
     &     -2.22327077943441828722D+01,
     &     +6.85986578867425578210D+00,
     &     -7.85596865167521363205D-01,
     &     +3.04288311387925105578D+01,
     &     -1.38697534163042433875D+02,
     &     +2.72764619820308155340D+02,
     &     -2.82250602598139892052D+02,
     &     +1.37077309352478550863D+02,
     &     -2.43363801583426564434D+01,
     &     -2.40391883810966504598D+02,
     &     +7.20123937632566253342D+02,
     &     -9.56047390568349555906D+02,
     &     +8.80126339460319968566D+02,
     &     -4.43060386994394036719D+02,
     &     +8.34534841118213392974D+01,
     &     +1.19040215556355087756D+03,
     &     -3.06320152050599153881D+03,
     &     +2.78453473219604893529D+03,
     &     -1.43582241909097660937D+03,
     &     +4.99720580375658244066D+02,
     &     -8.44179485778146130315D+01,
     &     -1.27279708339606577283D+03,
     &     +3.68563309227265335721D+03,
     &     -3.66200350857788362191D+03,
     &     +1.74828529177925724980D+03,
     &     -4.43446115665101842751D+02,
     &     +5.16038082522078980219D+01/
      
!     0.02000     LMA
      data (mrange(5,3,i),i=1,3) /0,9,1/
      data (rrange(5,3,i),i=1,3) /0,9,1/
      data (alphas(5,3,i),i=1,100) /
     &     -1.51049254914130983707D+03,
     &     +5.90804592923560267081D+03,
     &     -9.78411552163013584504D+03,
     &     +9.03281646038073631644D+03,
     &     -5.04733341717106668511D+03,
     &     +1.70840439977420237483D+03,
     &     -3.16251664756395882705D+02,
     &     +1.74202707579704139107D+01,
     &     +3.67713293499369164863D+00,
     &     -5.09055045574688058707D-01,
     &     +1.18642296961242536781D+05,
     &     -4.22882345998080738354D+05,
     &     +6.02728568485423107632D+05,
     &     -4.09618748060796875507D+05,
     &     +8.90979215331482701004D+04,
     &     +5.95414095173196401447D+04,
     &     -5.18590927506359803374D+04,
     &     +1.72331732723822569824D+04,
     &     -2.81277261133411593619D+03,
     &     +1.86487606226860776815D+02,
     &     -1.40854420295521779917D+06,
     &     +3.16920022432578727603D+06,
     &     +4.61300052596119523514D+05,
     &     -8.61105917229603044689D+06,
     &     +1.23556398056155946106D+07,
     &     -8.94036121833491511643D+06,
     &     +3.81112897666227677837D+06,
     &     -9.71158692551792715676D+05,
     &     +1.37401910653425060445D+05,
     &     -8.32876128728674302693D+03,
     &     +1.72458181532576158643D+07,
     &     -4.52772539129179418087D+07,
     &     +2.01148809934781342745D+07,
     &     +6.09227038679622411728D+07,
     &     -1.07679655401528432965D+08,
     &     +8.29801717072715312243D+07,
     &     -3.64493890278426855803D+07,
     &     +9.44902866139182634652D+06,
     &     -1.35153547941456316039D+06,
     &     +8.25342644314703647979D+04,
     &     -1.36692876826326847076D+08,
     &     +4.51025336598177015781D+08,
     &     -5.31011502233843207359D+08,
     &     +1.64806847028889596462D+08,
     &     +2.15493288928440541029D+08,
     &     -2.68977365806154668331D+08,
     &     +1.38776598747677832842D+08,
     &     -3.89913542043956518173D+07,
     &     +5.84939274544744472951D+06,
     &     -3.68464314099946233910D+05,
     &     +5.88007632627574801445D+08,
     &     -2.19546849066643857956D+09,
     &     +3.29578755621284484863D+09,
     &     -2.44640117334454774857D+09,
     &     +7.63562415504432678223D+08,
     &     +1.37457055430822163820D+08,
     &     -2.06529556624937295914D+08,
     &     +7.50013034409245699644D+07,
     &     -1.26735707259846013039D+07,
     &     +8.54653687692407402210D+05,
     &     -1.38327254361235642433D+09,
     &     +5.60614067713224506378D+09,
     &     -9.43342679709010314941D+09,
     &     +8.54102479910387897491D+09,
     &     -4.42476607837071800232D+09,
     &     +1.22015302616584563255D+09,
     &     -9.17636537403069436550D+07,
     &     -4.19402333332933112979D+07,
     &     +1.20023414749641995877D+07,
     &     -9.83676207032694481313D+05,
     &     +1.70433842457238912582D+09,
     &     -7.46564615829306602478D+09,
     &     +1.36213117744959812164D+10,
     &     -1.36449069424619369507D+10,
     &     +8.23382455328104782104D+09,
     &     -3.05178841976632213593D+09,
     &     +6.64362036529520750046D+08,
     &     -7.09507497840798497200D+07,
     &     +7.89097510117984027602D+05,
     &     +3.47019468185705947690D+05,
     &     -9.03133463246678471565D+08,
     &     +4.47621684993862724304D+09,
     &     -8.96349174100218391418D+09,
     &     +9.77159178772765159607D+09,
     &     -6.46222307786082077026D+09,
     &     +2.69448710911580371857D+09,
     &     -7.04153944276203036308D+08,
     &     +1.09216364710150659084D+08,
     &     -8.71251345965249091387D+06,
     &     +2.35094345904298679670D+05,
     &     +7.47826861721657812595D+07,
     &     -6.76538905537033677101D+08,
     &     +1.73090915363686656952D+09,
     &     -2.18473390758292007446D+09,
     &     +1.61084727759012627602D+09,
     &     -7.39637661531803250313D+08,
     &     +2.13872151698962450027D+08,
     &     -3.76329373731744587421D+07,
     &     +3.63438629714747751132D+06,
     &     -1.44132561512217595009D+05/
      
!     0.02000     HM
      data (mrange(5,4,i),i=1,3) /0,10,1/
      data (rrange(5,4,i),i=1,3) /-4,10,1/
      data (alphas(5,4,i),i=1,165) /
     &     +2.53729226092276048660D+09,
     &     -9.90007540851090812683D+09,
     &     +1.52593834113249797821D+10,
     &     -1.01549528384376564026D+10,
     &     +2.99829578010063916445D+07,
     &     +5.78893490124735832214D+09,
     &     -6.14609898490778446198D+09,
     &     +4.04278241615104389191D+09,
     &     -1.80654362590145325661D+09,
     &     +4.91334348953023314476D+08,
     &     -5.44061721887778788805D+07,
     &     -6.47383359300961066037D+06,
     &     +1.76558650155847985297D+06,
     &     +1.00096546886151510989D+05,
     &     -3.39961914804177431506D+04,
     &     -1.43765508369811935425D+10,
     &     +5.53386381611402206421D+10,
     &     -8.64972441405090179443D+10,
     &     +6.18263781388091888428D+10,
     &     -1.27031002818969459534D+10,
     &     -1.22741674009294738770D+10,
     &     +1.32545543793032627106D+10,
     &     -7.72103742705006504059D+09,
     &     +2.99797494889310359955D+09,
     &     -4.51316994491591095924D+08,
     &     -2.12359218253416419029D+08,
     &     +1.28106564253871873021D+08,
     &     -2.52927841362378038466D+07,
     &     +1.52832625190069107339D+06,
     &     +5.33991439763792805024D+04,
     &     +3.53888765741276016235D+10,
     &     -1.30949004957460281372D+11,
     &     +2.03740123625964263916D+11,
     &     -1.45485633356229827881D+11,
     &     +3.38169397075351562500D+10,
     &     +1.75117599864639778137D+10,
     &     -1.68526525631832427979D+10,
     &     +7.96671769498775196075D+09,
     &     -3.11421864218404531479D+09,
     &     +6.72407024087293863297D+08,
     &     +1.69479333713295638561D+08,
     &     -1.48711963495587915182D+08,
     &     +3.25634279867792949080D+07,
     &     -1.75731226257729995996D+06,
     &     -1.35079742770175478654D+05,
     &     -5.02259172277233886719D+10,
     &     +1.70340081645104858398D+11,
     &     -2.59178892024369842529D+11,
     &     +1.80456668127881774902D+11,
     &     -4.10109717020726623535D+10,
     &     -1.59603847360727577209D+10,
     &     +1.16582282107490043640D+10,
     &     -2.69128104888696050644D+09,
     &     +7.08683177906285405159D+08,
     &     -4.34911479088057935238D+08,
     &     +1.01044743218226939440D+08,
     &     +1.50179933915116582066D+07,
     &     -4.86460542015785723925D+06,
     &     -1.18190563420071476139D+06,
     &     +2.78623554417385661509D+05,
     &     +4.64705908425810546875D+10,
     &     -1.30875109087104644775D+11,
     &     +1.84384613263450225830D+11,
     &     -1.14315528459969711304D+11,
     &     +1.54155640604285240173D+10,
     &     +1.32641401056199913025D+10,
     &     -5.51530197341053485870D+09,
     &     -1.70800161550766736269D+08,
     &     +7.78830644995285868645D+08,
     &     -2.06842292132320821285D+08,
     &     -2.58832894926396012306D+07,
     &     +2.22092620407946482301D+07,
     &     -6.95661957240117993206D+06,
     &     +1.92690144883910729550D+06,
     &     -2.33277229786988784326D+05,
     &     -3.03660798334242973328D+10,
     &     +5.92434377201234664917D+10,
     &     -6.54677950676404571533D+10,
     &     +1.98006473040911102295D+10,
     &     +1.80741102193409042358D+10,
     &     -1.24110750642199649811D+10,
     &     +2.06936088729867076874D+09,
     &     -5.55200935775747060776D+08,
     &     +3.31793995932478666306D+08,
     &     -3.77282162409072965384D+07,
     &     +3.77271593664872169029D+05,
     &     -2.39515931142189726233D+06,
     &     +4.78976734628147096373D+05,
     &     -2.24492564827198133571D+05,
     &     +4.79792411042873136466D+04,
     &     +1.39020977034095058441D+10,
     &     -1.24214230056445674896D+10,
     &     +3.40159755029327201843D+09,
     &     +1.62527040079761524200D+10,
     &     -2.40572617307258148193D+10,
     &     +1.01576765639793376923D+10,
     &     -2.79329356560366511345D+08,
     &     -3.70225353787005186081D+08,
     &     -1.05583829363182082772D+08,
     &     +4.22791698590517193079D+07,
     &     -2.37510881991530815139D+06,
     &     -2.48551971484777750447D+06,
     &     +2.30980219498265953735D+06,
     &     -5.97702223240277264267D+05,
     &     +4.54666582402248095605D+04,
     &     -2.83252506215792465210D+09,
     &     -5.87107593619644260406D+09,
     &     +9.35640963237767410278D+09,
     &     -9.10900080282122421265D+09,
     &     +9.58587735224136924744D+09,
     &     -5.49030619645441913605D+09,
     &     +1.36275859178348875046D+09,
     &     -8.57005601597409099340D+07,
     &     -7.86904883546570241451D+07,
     &     +4.99692720418420732021D+07,
     &     -7.09569033777875266969D+06,
     &     -8.55338083130266284570D+05,
     &     -5.61264059584917849861D+05,
     &     +3.08073506627650989685D+05,
     &     -3.28288523737060459098D+04,
     &     -1.46525129680925726891D+09,
     &     +9.31328032341473007202D+09,
     &     -1.05576564874556827545D+10,
     &     +4.15304574768818807602D+09,
     &     -5.00505353118861734867D+08,
     &     +2.38979281271329969168D+08,
     &     -3.43075090044980287552D+08,
     &     +1.53399431073772042990D+08,
     &     +6.15683430398282781243D+06,
     &     -2.10229696774827726185D+07,
     &     +6.68991617223582346924D+05,
     &     +2.39900342987180408090D+06,
     &     -4.43480899599777883850D+05,
     &     -2.42983373542564295349D+04,
     &     +7.76394675149368777056D+03,
     &     +1.21311579654507946968D+09,
     &     -5.24070703990287113190D+09,
     &     +7.06700234139572620392D+09,
     &     -4.37710070889741706848D+09,
     &     +1.25442565132006978989D+09,
     &     -4.98846050957186818123D+07,
     &     -6.82854203030566722155D+07,
     &     +3.25620574902522787452D+07,
     &     -2.20218804932029694319D+07,
     &     +9.02162431089685671031D+06,
     &     -2.40492887758824275807D+05,
     &     -8.68585806973894243129D+05,
     &     +2.25402435254328796873D+05,
     &     -1.55503066541077005240D+04,
     &     -3.94672900362825259890D+02,
     &     -2.52924874342063546181D+08,
     &     +1.08433769070899581909D+09,
     &     -1.75274841405356764793D+09,
     &     +1.53857146715462470055D+09,
     &     -8.55304471302122354507D+08,
     &     +3.31495823209429860115D+08,
     &     -9.56754638305779546499D+07,
     &     +1.92680795062012895942D+07,
     &     -1.01626479920090991072D+06,
     &     -6.69255818740246933885D+05,
     &     +1.56196811545844484499D+04,
     &     +9.97457630392971332185D+04,
     &     -2.92981386761484645831D+04,
     &     +2.94914409981776043423D+03,
     &     -6.54138557526380992613D+01/
      
!     0.02000     Recom
      data (mrange(5,5,i),i=1,3) /0,3,1/
      data (rrange(5,5,i),i=1,3) /0,3,1/
      data (alphas(5,5,i),i=1,16) /
     &     +1.34142729153601880654D+01,
     &     -8.33977084929894418863D-01,
     &     +6.31803234438807592710D-01,
     &     -1.70217751115103232973D-01,
     &     +1.57297096125805579980D+00,
     &     +3.04879789211891349954D-01,
     &     -6.10228805816398045536D-01,
     &     +2.31090540934666771600D-01,
     &     -1.41108730414326188907D+00,
     &     +1.23004856793249794933D+00,
     &     -2.83350980192163703908D-01,
     &     -4.69140139332027000796D-02,
     &     +5.62788156340448986192D-01,
     &     -6.89182559641456582433D-01,
     &     +2.69502539954297459790D-01,
     &     -2.14046891159626988255D-02/