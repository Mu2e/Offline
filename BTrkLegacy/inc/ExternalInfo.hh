#ifndef BaBar_ExternalInfo
#define BaBar_ExternalInfo

//
// Manage access to externally supplied information.
//

class FileFinderInterface;
class ParticleInfoInterface;

class ExternalInfo {

public:

  // It is not possible to make an object of this type.
  // All information is accessed via static member functions.
  ExternalInfo  () = delete;
  ~ExternalInfo () = delete;
  ExternalInfo ( ExternalInfo const&  ) = delete;
  ExternalInfo ( ExternalInfo&&       ) = delete;
  ExternalInfo& operator= ( ExternalInfo const& ) = delete;
  ExternalInfo& operator= ( ExternalInfo&&      ) = delete;

  // Checked and unchecked accessors.
  static FileFinderInterface const* fileFinderInstance();
  static FileFinderInterface const* uncheckedFileFinderInstance(){
    return _findFile;
  }

  static ParticleInfoInterface const* particleInfoInstance();
  static ParticleInfoInterface const* uncheckedParticleInfoInstance(){
    return _particleInfo;
  }

  // Setters
  static void set ( FileFinderInterface const* findFile ){
    _findFile = findFile;
  }

  static void set ( ParticleInfoInterface const* particleInfo ){
    _particleInfo = particleInfo;
  }

private:

  // Non-owning pointers to externally managed objects.
  static FileFinderInterface   const* _findFile;
  static ParticleInfoInterface const* _particleInfo;

};

#endif // end BaBar_ExternalInfo
