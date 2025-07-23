#  Commande avancée d’un drone quadrirotor sous Matlab/Simulink

Projet réalisé dans le cadre du Master 1 – Université de Namur  

<p align="center">
  <img src="https://m.media-amazon.com/images/I/61Bd3hTy9jL.jpg" width="230" alt="MiniDrone"/>
  &nbsp;&nbsp;&nbsp;
  <img src="https://fiverr-res.cloudinary.com/images/q_auto,f_auto/gigs/354782021/original/5044f458d2b4544e3cfc9909a7d8c3ae051a5a25/do-matlab-coding-simulink-control-system-and-pid-controller-projects.png" width="320" alt="Simulink PID Control"/>
</p>


---

##  Description

Ce projet vise à modéliser, analyser et commander un drone quadrirotor à l’aide de Matlab et Simulink.  
Il comprend :
- **Étude complète du système** (modélisation dynamique, linéarisation, matrices d’état)
- **Analyse des propriétés** (stabilité, contrôlabilité, observabilité)
- **Synthèse et simulation de lois de commande** modernes :  
  - Placement de pôles  
  - Commande via LMIs (Linear Matrix Inequalities)  
  - Contrôleur PID
- **Application à un modèle réel** (Parrot MiniDrone Rolling Spider)

---

##  Fichiers du projet

- **Part2_Analyse** : code Matlab pour l’analyse dynamique du système
- **Part3_Asservissement.mat** : code Matlab pour tester le système en boucle ouverte/fermée (placement de pôles et LMI)
- **Feedback_control.slx** : fichier Simulink pour tester la loi de contrôle LMI sur le modèle non-linéaire
- **Controller_PID.slx** : fichier Simulink implémentant le contrôleur PID
- **Rapport_Sysco_Final.pdf** : rapport détaillé du projet (méthodo, résultats, analyse)

---

## Utilisation

1. Ouvrir les fichiers `.slx` sous Matlab/Simulink pour simuler les différentes stratégies de commande.
2. Utiliser les scripts Matlab (`Part2_Analyse`, `Part3_Asservissement.mat`) pour l’analyse et les tests.
3. Consulter le rapport PDF pour le détail théorique, les schémas et les résultats d’expérimentation.

---



##  Références

Voir `Rapport_Sysco_Final.pdf` pour la méthodologie complète, les équations, la validation et l’analyse critique.
