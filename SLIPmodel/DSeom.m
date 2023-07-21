function dx = DSeom(t,in1,in3,in4,Vl,Vs,h,Wi,l0,m,k,b,in13)
%DSeom
%    DX = DSeom(IN1,T,IN3,IN4,Vl,Vs,H,Wi,L0,M,K,B,IN13)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    03-Jul-2023 14:46:46

F1l = in3(:,1);
F2l = in3(:,2);
F1r = in4(:,1);
F2r = in4(:,2);
j_in1 = in13(:,1);
j_in2 = in13(:,2);
j_in3 = in13(:,3);
x1 = in1(:,1);
x2 = in1(:,2);
x3 = in1(:,3);
x4 = in1(:,4);
x5 = in1(:,5);
x6 = in1(:,6);
x7 = in1(:,7);
x8 = in1(:,8);
x9 = in1(:,9);
x10 = in1(:,10);
x11 = in1(:,11);
x12 = in1(:,12);
x13 = in1(:,13);
x14 = in1(:,14);
t2 = x1.*2.0;
t3 = x2.*2.0;
t4 = x3.*2.0;
t5 = x7.^2;
t6 = x8.^2;
t7 = x8.^3;
t8 = x9.^2;
t9 = x9.^3;
t10 = x10.^2;
t11 = x10.^3;
t12 = x11.^2;
t13 = x12.^2;
t14 = x13.^2;
t15 = x14.^2;
t16 = x7.*x8.*2.0;
t17 = x7.*x9.*2.0;
t18 = x7.*x10.*2.0;
t19 = x8.*x9.*2.0;
t20 = x8.*x10.*2.0;
t21 = x7.*x12.*2.0;
t22 = x8.*x11.*2.0;
t23 = x9.*x10.*2.0;
t24 = x7.*x13.*2.0;
t25 = x9.*x11.*2.0;
t26 = x7.*x12.*4.0;
t27 = x7.*x14.*2.0;
t28 = x8.*x11.*4.0;
t29 = x8.*x13.*2.0;
t30 = x9.*x12.*2.0;
t31 = x10.*x11.*2.0;
t32 = x7.*x13.*4.0;
t33 = x8.*x14.*2.0;
t34 = x9.*x11.*4.0;
t35 = x10.*x12.*2.0;
t36 = x7.*x14.*4.0;
t37 = x8.*x13.*4.0;
t38 = x9.*x12.*4.0;
t39 = x9.*x14.*2.0;
t40 = x10.*x11.*4.0;
t41 = x10.*x13.*2.0;
t42 = x8.*x14.*4.0;
t43 = x10.*x12.*4.0;
t44 = x9.*x14.*4.0;
t45 = x10.*x13.*4.0;
t46 = -F1l;
t47 = -F2l;
t48 = -F1r;
t49 = -j_in2;
t50 = -j_in3;
t51 = 1.0./j_in1;
t52 = 1.0./j_in2;
t53 = 1.0./j_in3;
t54 = 1.0./m;
t55 = -x1;
t57 = -x2;
t58 = -x3;
t60 = 1.0./x7;
t96 = j_in1.*x11.*x12.*8.0;
t97 = j_in1.*x11.*x13.*8.0;
t98 = j_in2.*x11.*x12.*8.0;
t99 = j_in1.*x11.*x14.*8.0;
t100 = j_in1.*x12.*x13.*8.0;
t101 = j_in2.*x11.*x13.*8.0;
t102 = j_in3.*x11.*x12.*8.0;
t103 = j_in1.*x12.*x14.*8.0;
t104 = j_in2.*x11.*x14.*8.0;
t105 = j_in2.*x12.*x13.*8.0;
t106 = j_in3.*x11.*x13.*8.0;
t107 = j_in1.*x13.*x14.*8.0;
t108 = j_in2.*x12.*x14.*8.0;
t109 = j_in3.*x11.*x14.*8.0;
t110 = j_in3.*x12.*x13.*8.0;
t111 = j_in2.*x13.*x14.*8.0;
t112 = j_in3.*x12.*x14.*8.0;
t113 = j_in3.*x13.*x14.*8.0;
t56 = -t2;
t59 = -t4;
t61 = t6.*2.0;
t62 = t8.*2.0;
t63 = t12.*2.0;
t64 = t13.*2.0;
t65 = t14.*2.0;
t66 = t15.*2.0;
t67 = -t19;
t68 = -t20;
t69 = -t22;
t70 = -t23;
t71 = -t25;
t72 = -t28;
t73 = -t30;
t74 = -t31;
t75 = -t33;
t76 = -t34;
t77 = -t37;
t78 = -t40;
t79 = -t41;
t80 = -t43;
t81 = -t44;
t82 = -t6;
t83 = -t8;
t84 = -t10;
t85 = F1l+t55;
t86 = F2l+t57;
t87 = F1r+t55;
t88 = F2r+t57;
t93 = j_in1+t49;
t94 = j_in1+t50;
t95 = j_in2+t50;
t114 = Vl.*t6.*-2.0;
t115 = Vl.*t8.*-2.0;
t116 = Vs.*t6.*-2.0;
t117 = Vs.*t8.*-2.0;
t118 = j_in1.*t12.*8.0;
t119 = j_in2.*t12.*8.0;
t120 = j_in1.*t14.*8.0;
t121 = j_in2.*t13.*8.0;
t122 = j_in3.*t12.*8.0;
t123 = j_in1.*t15.*8.0;
t124 = j_in3.*t13.*8.0;
t125 = j_in2.*t15.*8.0;
t126 = j_in3.*t14.*8.0;
t127 = -t98;
t128 = -t99;
t129 = -t100;
t130 = -t106;
t131 = -t111;
t132 = -t112;
t136 = t16+t23;
t137 = t17+t20;
t138 = t18+t19;
t89 = Vl.*t61;
t90 = Vl.*t62;
t91 = Vs.*t61;
t92 = Vs.*t62;
t133 = j_in1.*t95.*x8;
t134 = j_in2.*t94.*x9;
t135 = j_in3.*t93.*x10;
t139 = t16+t70;
t140 = t17+t68;
t141 = t18+t67;
t142 = j_in1.*t7.*t95;
t143 = j_in2.*t9.*t94;
t144 = j_in3.*t11.*t93;
t145 = Vl.*t137;
t146 = Vs.*t137;
t147 = h.*t137;
t149 = t49.*t94.*x9;
t150 = j_in1.*t6.*t95.*x7;
t151 = j_in2.*t8.*t94.*x7;
t152 = j_in3.*t10.*t93.*x7;
t157 = t61+t62-1.0;
t158 = (Wi.*t136)./2.0;
t168 = Vl+t114+t115+x3;
t169 = Vs+t116+t117+x3;
t174 = t5+t10+t82+t83;
t175 = t5+t8+t82+t84;
t176 = t5+t6+t83+t84;
t177 = t21+t39+t69+t79;
t178 = t24+t35+t71+t75;
t179 = t27+t29+t73+t74;
t180 = t26+t45+t72+t81;
t181 = t32+t42+t76+t80;
t182 = t36+t38+t77+t78;
t185 = t102+t107+t127;
t186 = t97+t108+t130;
t187 = t104+t110+t128;
t188 = t105+t109+t129;
t189 = t101+t103+t132;
t190 = t96+t113+t131;
t193 = t63+t64+t65+t66;
t195 = t120+t121+t122;
t196 = t119+t123+t124;
t197 = t118+t125+t126;
t148 = -t133;
t153 = Vl.*t139;
t154 = Vs.*t139;
t155 = h.*t139;
t156 = -t147;
t160 = t157.^2;
t161 = -t158;
t162 = (Wi.*t141)./2.0;
t164 = t46+t145+x1;
t165 = t48+t145+x1;
t166 = t46+t146+x1;
t167 = t48+t146+x1;
t183 = t168.^2;
t184 = t169.^2;
t191 = Vl.*t174;
t192 = Vs.*t174;
t194 = h.*t174;
t198 = t188.*x7;
t199 = t186.*x8;
t200 = t185.*x9;
t201 = t189.*x7;
t202 = t187.*x8;
t203 = t185.*x10;
t204 = t190.*x7;
t205 = t187.*x9;
t206 = t186.*x10;
t207 = j_in1.*t180.*x11;
t208 = j_in1.*t180.*x12;
t209 = j_in2.*t181.*x11;
t210 = j_in1.*t180.*x13;
t211 = j_in2.*t181.*x12;
t212 = j_in3.*t182.*x11;
t213 = j_in1.*t180.*x14;
t214 = j_in2.*t181.*x13;
t215 = j_in3.*t182.*x12;
t216 = j_in2.*t181.*x14;
t217 = j_in3.*t182.*x13;
t218 = j_in3.*t182.*x14;
t222 = t197.*x8;
t223 = t196.*x9;
t224 = t195.*x10;
t225 = (Wi.*t175)./2.0;
t230 = t50.*t182.*x12;
t231 = t49.*t181.*x14;
t159 = -t155;
t163 = -t162;
t170 = t86+t153;
t171 = t88+t153;
t172 = t86+t154;
t173 = t88+t154;
t219 = t191+x3;
t220 = t192+x3;
t221 = -t194;
t226 = -t201;
t227 = -t202;
t228 = -t203;
t229 = -t210;
t232 = -t222;
t233 = -t224;
t234 = t138.*t164;
t235 = t138.*t165;
t236 = t147+t162;
t239 = t156+t162;
t246 = t158+t194;
t247 = t166.*t175;
t248 = t167.*t175;
t249 = t161+t194;
t250 = t155+t225;
t267 = t160.*t183.*t184;
t275 = -t177.*(t155-t225);
t276 = -t179.*(t155-t225);
t279 = t208+t214+t218;
t283 = t209+t213+t230;
t284 = t207+t217+t231;
t237 = t141.*t172;
t238 = t141.*t173;
t242 = t85+t236;
t243 = t87+t147+t163;
t251 = t170.*t176;
t252 = t171.*t176;
t253 = t159+t225;
t254 = t158+t221+x3;
t255 = t58+t246;
t257 = t47+t250+x2;
t262 = t178.*t236;
t263 = t179.*t236;
t264 = -t178.*(t147+t163);
t265 = -t179.*(t147+t163);
t268 = t177.*t246;
t269 = t178.*t246;
t270 = -t177.*(t158+t221);
t271 = -t178.*(t158+t221);
t272 = t177.*t250;
t273 = t179.*t250;
t282 = t211+t212+t229;
t285 = t279.*x7.*2.0;
t287 = t284.*x8.*2.0;
t288 = t283.*x9.*2.0;
t240 = -t237;
t241 = -t238;
t244 = t242.^2;
t245 = t243.^2;
t256 = t254.^2;
t258 = t255.^2;
t259 = t88+t253;
t260 = t257.^2;
t266 = -t263;
t274 = -t269;
t277 = t234+t251;
t278 = t235+t252;
t286 = -t285;
t289 = t282.*x10.*2.0;
t292 = -t138.*t219.*(t237-t247);
t293 = -t138.*t219.*(t238-t248);
t294 = t4+t262+t272;
t295 = t3+t265+t268;
t301 = -t176.*t219.*(t237-t247);
t302 = -t176.*t219.*(t238-t248);
t303 = t59+t264+t275;
t304 = t56+t271+t273;
t261 = t259.^2;
t280 = t240+t247;
t281 = t241+t248;
t290 = t141.*t220.*t277;
t291 = t141.*t220.*t278;
t296 = t175.*t220.*t277;
t297 = t175.*t220.*t278;
t298 = t3+t266+t270;
t305 = t244+t256+t260;
t306 = t2+t274+t276;
t333 = t286+t287+t288+t289;
t299 = -t296;
t300 = -t297;
t307 = t245+t258+t261;
t308 = sqrt(t305);
t318 = t290+t301;
t319 = t291+t302;
t324 = (t296+t138.*t219.*(t237-t247)).^2;
t325 = (t297+t138.*t219.*(t238-t248)).^2;
t334 = t333.*x8;
t335 = t333.*x9;
t336 = t333.*x10;
t309 = 1.0./t308;
t310 = -t308;
t311 = sqrt(t307);
t320 = t292+t299;
t321 = t293+t300;
t322 = t318.^2;
t323 = t319.^2;
t337 = -t334;
t338 = -t336;
t312 = l0+t310;
t313 = 1.0./t311;
t315 = -t311;
t326 = t254.*t294.*t309;
t327 = t242.*t304.*t309;
t328 = -t257.*t309.*(-t3+t263+t177.*(t158+t221));
t339 = t267+t322+t324;
t340 = t267+t323+t325;
t314 = k.*t312;
t316 = l0+t315;
t329 = -t255.*t313.*(t4+t178.*(t147+t163)+t177.*(t155-t225));
t330 = t259.*t295.*t313;
t331 = -t243.*t313.*(t56+t269+t179.*(t155-t225));
t332 = t255.*t313.*(t4+t178.*(t147+t163)+t177.*(t155-t225));
t341 = 1.0./sqrt(t339);
t342 = 1.0./sqrt(t340);
t346 = t326+t327+t328;
t317 = k.*t316;
t343 = t157.*t168.*t169.*t254.*t309.*t341;
t344 = t157.*t168.*t169.*t255.*t313.*t342;
t347 = b.*t346;
t349 = t330+t331+t332;
t351 = t242.*t309.*t318.*t341;
t352 = t243.*t313.*t319.*t342;
t353 = -t257.*t309.*t341.*(t296+t138.*t219.*(t237-t247));
t354 = t257.*t309.*t341.*(t296+t138.*t219.*(t237-t247));
t355 = -t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248));
t345 = -t344;
t348 = -t347;
t350 = b.*t349;
t358 = t343+t351+t354;
t361 = -1.0./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248)));
t356 = t314+t348;
t357 = t317+t350;
t359 = 1.0./t358;
t360 = t345+t352+t355;
t362 = t85.*t157.*t168.*t169.*t341.*t356.*t359;
t363 = t86.*t157.*t168.*t169.*t341.*t356.*t359;
t366 = t87.*t157.*t168.*t169.*t342.*t357.*t361;
t367 = t88.*t157.*t168.*t169.*t342.*t357.*t361;
t368 = (t87.*t157.*t168.*t169.*t342.*t357)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248)));
t369 = (t88.*t157.*t168.*t169.*t342.*t357)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248)));
t370 = -t58.*t318.*t341.*t356.*t359;
t371 = t58.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247));
t372 = t86.*t318.*t341.*t356.*t359;
t373 = -t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247));
t375 = t319.*t342.*t357.*t361.*x3;
t376 = (t342.*t357.*x3.*(t297+t138.*t219.*(t238-t248)))./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248)));
t377 = t88.*t319.*t342.*t357.*t361;
t378 = (t87.*t342.*t357.*(t297+t138.*t219.*(t238-t248)))./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248)));
t379 = (t88.*t319.*t342.*t357)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248)));
t364 = -t362;
t365 = -t363;
t374 = -t372;
t382 = -t138.*(t362-t318.*t341.*t356.*t359.*x3);
t383 = -t139.*(t362-t318.*t341.*t356.*t359.*x3);
t384 = -t137.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)));
t385 = -t141.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)));
t386 = -t175.*(t362-t318.*t341.*t356.*t359.*x3);
t387 = -t176.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)));
t388 = t176.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)));
t389 = t368+t375;
t390 = t369+t376;
t391 = -t138.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))));
t392 = -t139.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))));
t395 = -t175.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))));
t399 = -t136.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)));
t400 = -t140.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)));
t401 = t136.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)));
t402 = -t174.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)));
t403 = t378+t379;
t380 = t364+t370;
t381 = t365+t371;
t393 = t137.*t390;
t394 = t141.*t390;
t396 = t176.*t390;
t398 = t373+t374;
t404 = t136.*t403;
t405 = t140.*t403;
t407 = t174.*t403;
t397 = -t396;
t406 = -t404;
t408 = t383+t384+t392+t393+t402+t407;
t409 = x7.*(-t393-t407+t139.*(t362-t318.*t341.*t356.*t359.*x3)+t139.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t174.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)))+t137.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)))).*-2.0;
t410 = x8.*(-t393-t407+t139.*(t362-t318.*t341.*t356.*t359.*x3)+t139.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t174.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)))+t137.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)))).*-2.0;
t411 = x9.*(-t393-t407+t139.*(t362-t318.*t341.*t356.*t359.*x3)+t139.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t174.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)))+t137.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)))).*-2.0;
t414 = x7.*(t387+t396-t405+t138.*(t362-t318.*t341.*t356.*t359.*x3)+t138.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t140.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)))).*-2.0;
t415 = x7.*(-t394+t399+t404+t175.*(t362-t318.*t341.*t356.*t359.*x3)+t175.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t141.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)))).*-2.0;
t416 = x8.*(-t394+t399+t404+t175.*(t362-t318.*t341.*t356.*t359.*x3)+t175.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t141.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)))).*-2.0;
t417 = x9.*(t387+t396-t405+t138.*(t362-t318.*t341.*t356.*t359.*x3)+t138.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t140.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)))).*-2.0;
t418 = x10.*(t387+t396-t405+t138.*(t362-t318.*t341.*t356.*t359.*x3)+t138.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t140.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)))).*-2.0;
t419 = x10.*(-t394+t399+t404+t175.*(t362-t318.*t341.*t356.*t359.*x3)+t175.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t141.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)))).*-2.0;
t420 = x7.*(t387+t396-t405+t138.*(t362-t318.*t341.*t356.*t359.*x3)+t138.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t140.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)))).*2.0;
t421 = x8.*(-t394+t399+t404+t175.*(t362-t318.*t341.*t356.*t359.*x3)+t175.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t141.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)))).*2.0;
t412 = t382+t388+t391+t397+t400+t405;
t413 = t385+t386+t394+t395+t401+t406;
t422 = t198+t199+t200+t233+t338+t409+t417+t421;
t423 = t204+t205+t206+t232+t337+t411+t419+t420;
t424 = t223+t226+t227+t228+t335+t410+t415+t418;
et1 = t201+t202+t203-t223-t335+x10.*(t387+t396-t405+t138.*(t362-t318.*t341.*t356.*t359.*x3)+t138.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t140.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)))).*2.0+x7.*(-t394+t399+t404+t175.*(t362-t318.*t341.*t356.*t359.*x3)+t175.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t141.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)))).*2.0;
et2 = x8.*(-t393-t407+t139.*(t362-t318.*t341.*t356.*t359.*x3)+t139.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t174.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)))+t137.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)))).*2.0;
et3 = t201+t202+t203-t223-t335+x10.*(t387+t396-t405+t138.*(t362-t318.*t341.*t356.*t359.*x3)+t138.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t140.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)))).*2.0+x7.*(-t394+t399+t404+t175.*(t362-t318.*t341.*t356.*t359.*x3)+t175.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t141.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)))).*2.0;
et4 = x8.*(-t393-t407+t139.*(t362-t318.*t341.*t356.*t359.*x3)+t139.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t174.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)))+t137.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)))).*2.0;
et5 = t201+t202+t203-t223-t335+x10.*(t387+t396-t405+t138.*(t362-t318.*t341.*t356.*t359.*x3)+t138.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t140.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)))).*2.0+x7.*(-t394+t399+t404+t175.*(t362-t318.*t341.*t356.*t359.*x3)+t175.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t141.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)))).*2.0;
et6 = x8.*(-t393-t407+t139.*(t362-t318.*t341.*t356.*t359.*x3)+t139.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t174.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)))+t137.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)))).*2.0;
et7 = t201+t202+t203-t223-t335+x10.*(t387+t396-t405+t138.*(t362-t318.*t341.*t356.*t359.*x3)+t138.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t140.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)))).*2.0+x7.*(-t394+t399+t404+t175.*(t362-t318.*t341.*t356.*t359.*x3)+t175.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t141.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)))).*2.0;
et8 = x8.*(-t393-t407+t139.*(t362-t318.*t341.*t356.*t359.*x3)+t139.*(t366+(t319.*t342.*t357.*x3)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))))+t174.*(t372+t85.*t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247)))+t137.*(t363+t341.*t356.*t359.*x3.*(t296+t138.*t219.*(t237-t247)))).*2.0;
mt1 = [x4;x5;x6;-t54.*(t318.*t341.*t356.*t359+t319.*t342.*t357.*t361);t54.*(t341.*t356.*t359.*(t296+t138.*t219.*(t237-t247))+t342.*t357.*t361.*(t297+t138.*t219.*(t238-t248)));-t54.*(m.*(9.81e+2./1.0e+2)-t157.*t168.*t169.*t341.*t356.*t359+(t157.*t168.*t169.*t342.*t357)./(t344-t352+t259.*t313.*t342.*(t297+t138.*t219.*(t238-t248))));x11;x12;x13;x14];
mt2 = [t193.*x7.*(-1.0./2.0)+(t51.*t52.*t53.*t60.*t423.*(j_in2.*j_in3.*x8+j_in3.*t8.*t93.*x8+j_in2.*t10.*t94.*x8+j_in1.*t95.*x7.*x9.*x10))./4.0+(t51.*t52.*t53.*t60.*t422.*(j_in1.*j_in2.*x10+j_in1.*t83.*t95.*x10+t6.*t49.*t94.*x10+j_in3.*t93.*x7.*x8.*x9))./4.0+(t51.*t52.*t53.*t60.*(j_in1.*j_in3.*x9+j_in1.*t10.*t95.*x9+t6.*t50.*t93.*x9+t49.*t94.*x7.*x8.*x10).*(et1+et2))./4.0;t193.*x8.*(-1.0./2.0)+(t51.*t52.*t53.*t60.*t422.*(t143+t149+t6.*t134+t135.*x7.*x8+j_in3.*t10.*t93.*x9))./4.0-(t51.*t52.*t53.*t60.*t423.*(t151+t152+j_in2.*j_in3.*x7+t133.*x9.*x10))./4.0-(t51.*t52.*t53.*t60.*(t144+t6.*t135+t50.*t93.*x10+t149.*x7.*x8+j_in2.*t8.*t94.*x10).*(et3+et4))./4.0];
mt3 = [t193.*x9.*(-1.0./2.0)-(t51.*t52.*t53.*t60.*t422.*(t142+t148+t8.*t133+t135.*x7.*x9+t10.*t50.*t93.*x8))./4.0-(t51.*t52.*t53.*t60.*(t150+j_in1.*j_in3.*x7+t149.*x8.*x10+t10.*t50.*t93.*x7).*(et5+et6))./4.0+(t51.*t52.*t53.*t60.*t423.*(t135+t11.*t50.*t93+t133.*x7.*x9+j_in1.*t6.*t95.*x10+t8.*t50.*t93.*x10))./4.0;t193.*x10.*(-1.0./2.0)+(t51.*t52.*t53.*t60.*t423.*(t143+t149+t10.*t134+t148.*x7.*x10+j_in1.*t6.*t95.*x9))./4.0-(t51.*t52.*t53.*t60.*(t142+t148+t10.*t133+t134.*x7.*x10+j_in2.*t8.*t94.*x8).*(et7+et8))./4.0-(t51.*t52.*t53.*t60.*t422.*(j_in1.*j_in2.*x7+t135.*x8.*x9+j_in1.*t82.*t95.*x7+t8.*t49.*t94.*x7))./4.0];
dx = [mt1;mt2;mt3];
