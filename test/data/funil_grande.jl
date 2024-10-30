# Funil grande from 1931 to 2019
const funil_grande = [302.,537.,298.,167.,80.0,43.4,85.0,76.0,75.0,35.0,76.0,158.,
159.,261.,282.,220.,163.,105.,84.0,71.0,73.0,96.0,110.,327.,
407.,198.,178.,119.,86.0,72.0,68.0,62.0,59.0,90.0,91.0,169.,
278.,129.,141.,86.0,67.0,58.0,53.0,46.0,51.0,57.0,81.0,252.,
383.,539.,263.,227.,158.,115.,106.,81.0,69.0,100.0,100.0,117.,
81.0,102.,259.,163.,92.0,67.0,58.0,59.0,62.0,58.0,129.,186.,
540.,333.,227.,155.,137.,98.0,78.0,63.0,57.0,163.,253.,561.,
367.,397.,329.,227.,182.,129.,100.0,91.0,89.0,144.,199.,404.,
466.,393.,231.,207.,138.,101.,85.0,67.0,60.0,81.0,82.0,256.,
321.,432.,340.,214.,139.,102.,81.0,65.0,66.0,79.0,217.,269.,
444.,270.,267.,208.,135.,116.,110.,85.0,107.,124.,158.,344.,
327.,263.,372.,207.,160.,122.,101.,86.0,80.0,118.,182.,341.,
734.,556.,511.,320.,204.,111.,96.0,79.0,74.0,96.0,91.0,216.,
194.,342.,327.,180.,119.,84.0,72.0,58.0,50.0,55.0,71.0,128.,
184.,286.,191.,177.,102.,83.0,68.0,53.0,51.0,55.0,93.0,243.,
460.,180.,209.,164.,110.,87.0,73.0,63.0,56.0,82.0,148.,143.,
270.,261.,686.,286.,186.,142.,115.,100.0,119.,101.,108.,268.,
270.,280.,298.,181.,125.,99.0,79.0,65.0,58.0,57.0,111.,281.,
298.,399.,237.,167.,113.,120.,83.0,64.0,56.0,75.0,76.0,231.,
269.,392.,240.,189.,141.,106.,86.0,70.0,65.0,83.0,261.,327.,
294.,447.,396.,273.,168.,135.,110.,94.0,79.0,88.0,74.0,160.,
274.,369.,415.,259.,158.,136.,106.,85.0,71.0,76.0,138.,212.,
152.,149.,145.,185.,113.,89.0,69.0,46.0,61.0,61.0,98.0,170.,
109.,173.,101.,125.,84.0,45.0,54.0,46.0,41.0,54.0,90.0,127.,
172.,129.,137.,110.,71.0,63.0,40.0,33.1,27.7,33.0,47.0,158.,
188.,135.,148.,107.,88.0,78.0,60.0,61.0,53.0,33.0,60.0,244.,
218.,193.,213.,249.,137.,117.,93.0,74.0,72.0,54.0,97.0,220.,
117.,115.,153.,115.,121.,82.0,66.0,50.0,56.0,44.0,40.0,90.0,
257.,201.,142.,146.,126.,112.,98.0,79.0,70.0,84.0,100.0,131.,
285.,321.,404.,177.,158.,120.,135.,106.,76.0,81.0,95.0,225.,
453.,497.,480.,265.,234.,160.,132.,114.,106.,85.0,125.,144.,
256.,361.,311.,181.,161.,153.,117.,103.,92.0,146.,159.,231.,
329.,327.,188.,124.,97.0,78.0,71.0,59.0,52.0,73.0,100.0,74.1,
316.,392.,145.,172.,135.,103.,96.0,80.0,68.0,105.,145.,246.,
509.,671.,499.,336.,286.,217.,192.,158.,179.,258.,313.,344.,
607.,310.,362.,228.,148.,119.,103.,88.0,79.0,185.,270.,367.,
550.,416.,285.,186.,146.,130.,106.,97.0,78.0,92.0,231.,244.,
298.,189.,191.,124.,98.0,81.0,74.0,78.0,90.0,111.,105.,311.,
270.,255.,208.,148.,98.0,108.,89.0,81.0,74.0,135.,272.,259.,
300.,207.,172.,144.,95.0,79.0,80.0,74.0,85.0,121.,191.,131.,
150.,96.0,102.,74.0,65.0,82.0,67.0,63.0,67.0,96.0,168.,373.,
211.,290.,310.,174.,129.,102.,118.,92.0,86.0,126.,300.,280.,
341.,332.,155.,195.,148.,110.,106.,100.0,97.0,115.,170.,257.,
320.,228.,283.,238.,147.,111.,96.0,82.0,66.0,102.,102.,208.,
294.,265.,154.,115.,98.0,81.0,84.0,70.0,60.0,98.0,214.,257.,
174.,184.,205.,138.,122.,110.,114.,118.,169.,171.,226.,330.,
283.,240.,223.,192.,121.,109.,94.0,80.0,103.,85.0,160.,251.,
495.,226.,200.,156.,127.,120.,104.,86.0,85.0,110.,187.,221.,
323.,684.,329.,184.,160.,134.,127.,109.,137.,120.,190.,338.,
509.,321.,191.,238.,144.,126.,111.,96.0,89.0,95.0,152.,351.,
420.,229.,223.,144.,119.,113.,88.0,92.0,81.0,120.,276.,352.,
402.,259.,424.,268.,167.,142.,123.,113.,91.0,160.,206.,425.,
625.,462.,508.,377.,222.,217.,182.,151.,181.,283.,341.,551.,
318.,194.,180.,145.,140.,102.,81.0,71.0,96.0,94.0,129.,303.,
515.,422.,434.,219.,155.,128.,108.,90.0,95.0,103.,163.,221.,
354.,305.,225.,139.,122.,92.0,90.0,100.0,79.0,62.0,74.0,353.,
349.,310.,231.,175.,143.,120.,99.0,75.0,97.0,88.0,112.,315.,
301.,416.,249.,168.,133.,114.,87.0,75.0,63.0,116.,160.,206.,
341.,305.,308.,174.,123.,116.,97.0,93.0,94.0,75.0,89.0,205.,
251.,128.,199.,146.,119.,91.0,89.0,79.0,97.0,94.0,111.,134.,
485.,367.,296.,260.,153.,121.,100.0,81.0,72.0,118.,123.,172.,
803.,467.,231.,171.,145.,100.0,82.0,73.0,95.0,110.,231.,182.,
234.,302.,289.,217.,133.,126.,92.0,78.0,77.0,123.,103.,159.,
451.,208.,282.,183.,191.,130.,103.,83.0,69.0,79.0,95.0,216.,
172.,400.,215.,161.,129.,103.,87.0,65.0,63.0,91.0,130.,215.,
349.,249.,241.,162.,122.,106.,95.0,74.0,99.0,96.0,312.,360.,
839.,318.,305.,188.,137.,130.,104.,84.0,77.0,85.0,123.,210.,
230.,235.,179.,133.,105.,96.0,75.0,76.0,59.0,86.0,133.,212.,
264.,211.,279.,151.,105.,91.0,68.0,54.0,45.0,49.0,78.0,162.,
329.,305.,224.,152.,101.,89.0,77.0,66.0,90.0,63.0,129.,171.,
172.,125.,123.,116.,75.0,61.0,52.0,49.0,54.0,67.0,129.,241.,
246.,348.,205.,116.,87.0,96.0,84.0,67.0,62.0,43.0,87.0,134.,
329.,223.,193.,114.,88.0,67.0,58.0,54.0,51.0,45.0,75.0,198.,
263.,296.,254.,208.,123.,110.,94.0,70.0,73.0,62.0,84.0,307.,
428.,273.,278.,161.,145.,115.,93.0,74.0,75.0,64.0,149.,341.,
190.,222.,220.,124.,102.,83.0,71.0,58.0,64.0,96.0,130.,215.,
498.,307.,162.,124.,98.0,86.0,73.0,56.0,43.0,54.0,92.0,160.,
197.,308.,339.,232.,128.,103.,82.0,71.0,81.0,90.0,159.,391.,
454.,355.,299.,304.,152.,127.,104.,79.0,97.0,160.,151.,351.,
345.,191.,298.,172.,124.,102.,88.0,67.0,61.0,98.0,211.,326.,
430.,165.,321.,197.,120.,112.,87.0,68.0,50.0,79.0,100.0,417.,
685.,274.,222.,188.,149.,138.,108.,85.0,69.0,68.0,102.,126.,
270.,250.,207.,178.,113.,111.,83.0,62.0,63.0,85.0,95.0,213.,
123.,72.0,79.0,86.0,56.0,50.0,44.0,41.0,31.0,24.0,88.0,111.,
75.4,136.,187.,118.,83.0,64.0,50.0,40.0,75.0,38.0,118.,201.,
288.,210.,187.,96.0,78.0,82.0,58.0,48.0,46.0,54.0,141.,174.,
155.,115.,106.,74.0,74.0,63.0,47.0,37.0,29.0,46.0,60.0,144.,
170.,141.,210.,88.0,59.0,55.0,44.0,57.0,43.0,73.0,135.,197.,
134.,145.,215.,127.,87.0,69.0,54.0,47.0,40.0,45.0,100.0,158.]