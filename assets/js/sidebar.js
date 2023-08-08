const addedTopBounding = 80;
const headerOffSetToAdd = 7;
document.addEventListener("DOMContentLoaded", function () {
  copySidebar()
  addTopMobileNav()
  const sidebarToggle = document.querySelector(".sidebar-nav");
  const sidebarContentLinks = document.querySelectorAll(".internal-nav a");
  const headings = document.querySelectorAll("h1, h2, h3, h4, h5, h6");
  const header = document.querySelector(".header-size");
  const headerHeight = header.offsetHeight + headerOffSetToAdd


  sidebarToggle?.addEventListener("click", () => {
    toggleNavMenu(sidebarToggle);
  });

  sidebarContentLinks.forEach((link) => {
    link.addEventListener("click", function (event) {
      event.preventDefault();

      const target = document.querySelector(this.getAttribute("href"));
      if (target) {
        const offsetTop = target.offsetTop - headerHeight
        window.scrollTo({
          top: offsetTop,
          behavior: "smooth",
        });
      }
    });
  });

  window.addEventListener("scroll", () => {
    setCurrentHeading(sidebarContentLinks, findCurrentHeadingId(headings));
  });
});


function findCurrentHeadingId(headings) {
  for (const heading of headings) {
    const bounding = heading.getBoundingClientRect();
    if (heading.id) {
      if (
        bounding.top >= addedTopBounding &&
        bounding.bottom <= window.innerHeight
      ) {
        return heading.id;
      }
    }
  }
}

function setCurrentHeading(sidebarContentlinks, activeHeading) {
  sidebarContentlinks.forEach((link) => {
    const linkHeading = link.getAttribute("href").slice(1);
    activeHeading === linkHeading
      ? link.classList.add("selected-nav")
      : link.classList.remove("selected-nav");
  });
}

function toggleNavMenu(sidebarContent) {
  sidebarContent.classList.contains("open")
    ? sidebarContent.classList.remove("open")
    : sidebarContent.classList.add("open");
}

function copySidebar() {

  const header = document.querySelector("header")
  const sidebarContainer = document.querySelector(".sidebar-nav-container");
  
  if(sidebarContainer){
    const sidebarCopy = sidebarContainer.cloneNode(true);
    sidebarCopy.id = "header-sidebar-container";

    header.appendChild(sidebarCopy);
  }
  
}

function addTopMobileNav() {
  
  const sideBarNav = document.querySelector("#header-sidebar-container");
  const headerNav = document.querySelector(".header-nav");
  sideBarNav ?
    headerNav.style.top = "60%" : headerNav.style.top = "100%"
   
}
